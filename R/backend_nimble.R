#' Compile and run a NIMBLE model with standard because settings
#'
#' @param nimble_code NIMBLE model code (from nimbleCode)
#' @param constants Data constants for NIMBLE
#' @param data NIMBLE data list
#' @param inits Initial values list
#' @param n.chains Number of MCMC chains
#' @param n.iter Number of iterations
#' @param n.burnin Number of burnin iterations
#' @param n.thin Thinning interval
#' @param monitor Character vector of parameters to monitor
#' @param quiet Logical, suppress output
#' @param family Named vector of family types for hardening samplers
#' @param nimble_samplers Optional custom sampler list
#' @param nimble_waic Logical, calculate WAIC
#' @return A list containing samples and the compiled model
#' @export
run_nimble_model <- function(nimble_code, constants, data, inits, 
                             n.chains, n.iter, n.burnin, n.thin, 
                             monitor, quiet, family = NULL, 
                             nimble_samplers = NULL, nimble_waic = TRUE) {
  
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop("The 'nimble' package is required to use engine = 'nimble'")
  }
  
  model <- tryCatch({
    nimble::nimbleModel(code = nimble_code, constants = constants, data = data, inits = inits[[1]])
  }, error = function(e) {
    if (!quiet) {
      message("\nCRITICAL NIMBLE ERROR during model initialization:")
      message(e$message)
    }
    stop(e)
  })
  
  cModel <- tryCatch({
    nimble::compileNimble(model, showCompilerOutput = !quiet)
  }, error = function(e) {
    stop(paste("NIMBLE compilation failed:", e$message))
  })
  
  mcmc_conf <- nimble::configureMCMC(model, monitors = monitor, print = !quiet)
  
  # Apply because sampler hardening
  if (exists("nimble_harden_samplers", mode = "function")) {
    nimble_harden_samplers(mcmc_conf, family, nimble_samplers, quiet)
  }
  
  mcmc <- nimble::buildMCMC(mcmc_conf)
  cMcmc <- nimble::compileNimble(mcmc, project = model, showCompilerOutput = !quiet)
  
  samples <- nimble::runMCMC(
    cMcmc, 
    niter = n.iter, 
    nburnin = n.burnin, 
    nchains = n.chains, 
    thin = n.thin, 
    inits = inits,
    samplesAsCodaMCMC = TRUE, 
    WAIC = nimble_waic,
    summary = FALSE
  )
  
  # Format output to match standard because structure
  result <- list()
  if (n.chains == 1) {
    result$samples <- coda::mcmc.list(coda::as.mcmc(samples))
  } else {
    if (nimble_waic) {
      result$samples <- coda::mcmc.list(lapply(samples$samples, coda::as.mcmc))
      result$WAIC <- samples$WAIC
    } else {
      result$samples <- coda::mcmc.list(lapply(samples, coda::as.mcmc))
    }
  }
  
  result$model <- cModel
  return(result)
}

#' Harden NIMBLE Sampler Configuration
#'
#' Applies robust sampler assignments to a NIMBLE MCMC configuration object.
#' Exported as a standalone function so that parallel worker nodes --- which
#' load the package fresh --- always use the current installed version of this
#' logic, regardless of which version of `because()` was originally called.
#'
#' @param mcmc_conf A NIMBLE MCMC configuration object (from `configureMCMC()`).
#' @param family Named character vector of response families (same as `because()`).
#' @param nimble_samplers Optional named list of user-specified samplers.
#' @param quiet Logical. If TRUE, suppress status messages.
#' @return The modified `mcmc_conf` object (invisibly).
#' @keywords internal
#' @export
nimble_harden_samplers <- function(mcmc_conf, family = NULL, nimble_samplers = NULL, quiet = TRUE) {

  # Helper: parse VAR name from an err_raw_VAR_STRUCTURE[...] or u_std_VAR_STRUCTURE[...] node
  # Strategy: strip prefix and index, then the last _-delimited token is the structure name.
  parse_re_info <- function(node) {
    base   <- sub("\\[.*\\]", "", node)                  # strip [1:N]
    base   <- sub("^(err_raw_|u_std_|sigma_|tau_|beta_|alpha_)", "", base) # strip prefixes
    parts  <- strsplit(base, "_")[[1]]
    if (length(parts) < 2) return(list(var = base, structure = ""))
    
    # For betas: beta_RESPONSE_PREDICTOR
    if (grepl("^beta_", node)) {
       return(list(var = parts[1], structure = "beta"))
    }

    # For REs and scales: name_STRUCTURE
    last   <- parts[length(parts)]
    known_structs <- c("phylo", "spatial", "survey", "site", "obs", "res")
    if (last %in% known_structs) {
       return(list(
         var       = paste(parts[-length(parts)], collapse = "_"), 
         structure = last                                          
       ))
    }
    list(var = base, structure = "")
  }

  sampler_targets <- sapply(mcmc_conf$getSamplers(), function(x) x$target)

  # ------ 1. Hybrid Grouping Logic: Core Equation Blocks ------------------------------------------------------------------------
  # The goal is a SMALL block (5-10 nodes) containing: (alpha, all betas, all scales)
  trait_groups <- list()

  # Identify possible intercepts
  alpha_nodes <- unique(grep("^alpha_.*", sampler_targets, value = TRUE))
  for (a_node in alpha_nodes) {
    trait_name <- sub("^alpha_", "", a_node)
    trait_groups[[trait_name]] <- list(targets = a_node)
  }

  # Identify all Fixed Effects (Slopes) and Variance nodes
  hyper_nodes <- unique(grep("^(beta_|sigma_|tau_).*", sampler_targets, value = TRUE))
  
  for (node in hyper_nodes) {
    # Skip if already handled by posterior_predictive logic
    current_types <- sapply(mcmc_conf$getSamplers(node), function(s) s$name)
    if (any(grepl("posterior_predictive", current_types, ignore.case = TRUE))) next

    info <- parse_re_info(node)
    t_name <- info$var
    
    # Heuristic for beta_RESPONSE_PREDICTOR: find the matching trait name
    if (grepl("^beta_", node)) {
       # Find which trait_group it belongs to
       matches <- names(trait_groups)[vapply(names(trait_groups), function(tn) grepl(paste0("^", tn, "_"), sub("^beta_", "", node)), logical(1))]
       if (length(matches) > 0) t_name <- matches[which.max(nchar(matches))]
    }

    # Clean up standard trait name prefixing
    if (t_name == "" || is.na(t_name)) next
    if (grepl("^phylo_", t_name)) t_name <- sub("^phylo_", "", t_name)

    # Add to group if a matching intercept (alpha) exists
    if (!is.null(trait_groups[[t_name]])) {
       trait_groups[[t_name]]$targets <- unique(c(trait_groups[[t_name]]$targets, node))
    }
  }

  # ------ 2. Apply Joint Core Blocks (AF_slice vs RW_block) ---------------------------------------------------------------------------
  processed_nodes <- character(0)

  for (t_name in names(trait_groups)) {
    group_targets <- trait_groups[[t_name]]$targets
    
    # We block if there's an intercept AND at least one scale/slope component
    if (length(group_targets) > 1) {
      for (target in group_targets) {
        mcmc_conf$removeSamplers(target)
      }

      # Root Node Detection: Traits with no slopes (betas) get AF_slice
      has_slopes   <- any(grepl("^beta_", group_targets))
      sampler_type <- if (has_slopes) "RW_block" else "AF_slice"

      mcmc_conf$addSampler(
        target = group_targets,
        type   = sampler_type
      )
      
      processed_nodes <- c(processed_nodes, group_targets)

      if (!quiet) {
        message(sprintf(
          "NIMBLE: %s Core Block for '%s' (%d nodes: %s, etc.)",
          sampler_type, t_name, length(group_targets), group_targets[1]
        ))
      }
    }
  }

  # ------ 3. Handle Remaining Nodes (ESS and Defaults) ------------------------------------------------------------------------------
  remaining_targets <- setdiff(sampler_targets, processed_nodes)
  
  for (target in remaining_targets) {
    current_types <- sapply(mcmc_conf$getSamplers(target), function(s) s$name)
    if (any(grepl("posterior_predictive", current_types, ignore.case = TRUE))) next
    
    # [NEW] Phylogenetic ESS: Force Elliptical Slice Sampler for multivariate phylo nodes.
    # Matches both centered (err_raw_*_phylo) and non-centered (z_*_phylo) parameterizations.
    if (grepl("^(err_raw_|u_std_|z_).*(phylo|BM|OU|Pagel).*", target)) {
        mcmc_conf$removeSamplers(target)
        mcmc_conf$addSampler(target = target, type = "ess")
        if (!quiet) message(sprintf("NIMBLE: Forced 'ess' (Elliptical Slice) for phylogenetic vector '%s'.", target))
        next
    }

    # Other REs (Standalone / Default)
    if (grepl("^(err_raw_|u_std_).*", target)) {
        if (length(current_types) > 1 || (!any(grepl("RW_block|conjugate|ess", current_types, ignore.case = TRUE)))) {
            mcmc_conf$removeSamplers(target)
            mcmc_conf$addSampler(target = target, type = "RW_block")
        }
    }
    
    # Standalone Scales (Standalone Defaults)
    if (grepl("^(sigma_|tau_|lambda_|r_|psi_|sigma_total_).*", target)) {
      mcmc_conf$removeSamplers(target)
      mcmc_conf$addSampler(target = target, type = "slice")
    }
  }

  if (!quiet) {
      message(sprintf("NIMBLE: Sampler hardening complete. ESS + AF_slice + Hybrid strategy applied."))
  }

  # ------ 7. User-specified overrides (always applied last) ---------------------------------------------------------------
  if (!is.null(nimble_samplers)) {
    for (node in names(nimble_samplers)) {
      mcmc_conf$removeSamplers(node)
      mcmc_conf$addSampler(target = node, type = nimble_samplers[[node]])
    }
  }

  invisible(mcmc_conf)
}

#' Execute full NIMBLE pipeline (sequential or parallel) for because()
#'
#' @keywords internal
run_nimble_pipeline <- function(
  model_string, data, family, extension_inits, monitor,
  n.chains, n.iter, n.burnin, n.thin, WAIC,
  nimble_samplers, parallel, n.cores, cl = NULL, quiet = FALSE
) {
  if (!requireNamespace("nimble", quietly = TRUE)) {
    stop(
      "The 'nimble' package is required when engine = 'nimble'.\n",
      "Please install it using: install.packages('nimble')\n",
      "For detailed installation instructions and system requirements (e.g. Rtools/Xcode),\n",
      "see: https://r-nimble.org/download"
    )
  }

  if (!("package:nimble" %in% search())) {
    suppressPackageStartupMessages(attachNamespace("nimble"))
  }

  if (!quiet) {
    message("Compiling model via NIMBLE...")
  }

  nimble_funcs <- list()
  nimble_string <- model_string

  # Generic cleanup for NIMBLE: strip JAGS-specific log-density nodes
  lines <- strsplit(nimble_string, "\n")[[1]]
  lines <- lines[!grepl("logdensity\\.", lines)]
  lines <- lines[!grepl("log_lik_", lines)]
  lines <- lines[!grepl("lik_matrix_", lines)]
  nimble_string <- paste(lines, collapse = "\n")

  for (v in names(family)) {
    fam_obj <- structure(
      list(name = family[[v]]),
      class = c(paste0("because_family_", family[[v]]), "because_family")
    )
    opt_res <- nimble_family_optimization(
      fam_obj,
      nimble_string,
      variable = v
    )
    nimble_string <- opt_res$model_string
    if (length(opt_res$nimble_functions) > 0) {
      nimble_funcs <- c(nimble_funcs, opt_res$nimble_functions)
    }

    # If discretized latent states were marginalized, remove them from monitors
    if (
      nimble_string != model_string && any(grepl(paste0("z_", v), monitor))
    ) {
      monitor <- setdiff(monitor, paste0("z_", v))
    }
  }

  # Register nimble functions to local environment for compiler
  if (length(nimble_funcs) > 0) {
    unique_names <- unique(names(nimble_funcs))
    for (fn_name in unique_names) {
      assign(fn_name, nimble_funcs[[fn_name]], envir = environment())
      if (startsWith(fn_name, "d")) {
        try(
          nimble::registerDistributions(nimble_funcs[fn_name]),
          silent = TRUE
        )
      }
    }
  }

  nimble_inits <- extension_inits
  if (is.null(nimble_inits)) nimble_inits <- list()
  for (p in monitor) {
    if (!p %in% names(nimble_inits)) {
      if (grepl("^(tau_|sigmay_|sigmap_|sigmar_|sigma_)", p)) {
        nimble_inits[[p]] <- 1.0
      } else if (grepl("^alpha_", p)) {
        resp_name <- sub("^alpha_", "", p)
        if (resp_name %in% names(data)) {
            m_val <- mean(as.numeric(data[[resp_name]]), na.rm = TRUE)
            if (all(as.numeric(data[[resp_name]]) >= 0, na.rm = TRUE)) {
                nimble_inits[[p]] <- log(max(0.1, m_val))
            } else {
                nimble_inits[[p]] <- m_val
            }
        } else {
          nimble_inits[[p]] <- 0.0
        }
      } else if (grepl("^beta_", p)) {
        nimble_inits[[p]] <- 0.0
      } else if (grepl("^psi_", p)) {
        nimble_inits[[p]] <- 0.5
      } else if (grepl("^r_", p)) {
        nimble_inits[[p]] <- 1.0
      } else if (grepl("^sigma_total_", p)) {
        nimble_inits[[p]] <- 1.0
      } else if (grepl("^lambda_", p)) {
        nimble_inits[[p]] <- 0.5
      } else if (grepl("^cutpoint", p)) {
        nimble_inits[[p]] <- 0.0
      }
    }
  }

  nimble_string_clean <- sub("^\\s*model\\s*\\{", "{", nimble_string)
  nimble_code <- parse(text = nimble_string_clean)[[1]]

  nimble_constants <- data
  nimble_data <- list()
  if (!is.null(data[["L_multiPhylo"]])) {
      nimble_data[["L_multiPhylo"]] <- data[["L_multiPhylo"]]
      nimble_constants[["L_multiPhylo"]] <- NULL
  }
  if (!is.null(data[["Prec_multiPhylo"]])) {
      nimble_data[["Prec_multiPhylo"]] <- data[["Prec_multiPhylo"]]
      nimble_constants[["Prec_multiPhylo"]] <- NULL
  }
  if (!is.null(data[["L_phylo"]])) {
      nimble_data[["L_phylo"]] <- data[["L_phylo"]]
      nimble_constants[["L_phylo"]] <- NULL
  }

  nimble_model <- tryCatch(
    {
      m_obj <- suppressMessages(suppressWarnings(nimble::nimbleModel(
        code = nimble_code,
        constants = nimble_constants,
        data = nimble_data,
        inits = nimble_inits,
        buildDerivs = (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC"),
        calculate = FALSE
      )))
      
      model_nodes <- m_obj$getNodeNames(stochOnly = TRUE, includeData = FALSE)
      for (node in model_nodes) {
        base_node <- sub("\\[.*\\]", "", node)
        if (!any(grepl(paste0("^", base_node, "$"), names(nimble_inits)))) {
          is_prec <- grepl("^(tau_|sigmay_|sigmap_|sigmar_)", node)
          val <- if (is_prec) 1 else 0
          if (grepl("^psi_", node)) val <- 0.5
          if (grepl("^r_", node)) val <- 1
          if (grepl("^sigma_total_", node)) val <- 1.0
          if (grepl("^lambda_", node)) val <- 0.5
          
          try({
            curr_val <- m_obj[[node]]
            if (any(is.na(curr_val)) || any(is.nan(curr_val))) {
              m_obj[[node]] <- val
            }
          }, silent = TRUE)
        }
      }
      
      unique_base_nodes <- unique(sub("\\[.*\\]", "", model_nodes))
      for (v in unique_base_nodes) {
          if (!v %in% names(nimble_inits)) {
              try({
                  nimble_inits[[v]] <- m_obj[[v]]
              }, silent = TRUE)
          }
      }
      
      m_obj
    },
    error = function(e) {
      if (!quiet) {
        message("\nCRITICAL NIMBLE ERROR during model initialization:")
        message(e$message)
      }
      stop(paste(e, "\n\n", model_string))
    }
  )

  mcmc_conf <- nimble::configureMCMC(
    nimble_model,
    monitors = monitor,
    enableWAIC = WAIC
  )

  if (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC") {
    if (!requireNamespace("nimbleHMC", quietly = TRUE)) {
      stop("The 'nimbleHMC' package is required to use HMC samplers in NIMBLE. Install it with install.packages('nimbleHMC')")
    }
    if (!("package:nimbleHMC" %in% search())) {
      suppressPackageStartupMessages(attachNamespace("nimbleHMC"))
    }
    if (!quiet) {
      message("Building and compiling NIMBLE MCMC using nimbleHMC::buildHMC (this may take a moment)...")
    }
    nimble_mcmc <- nimbleHMC::buildHMC(nimble_model)
  } else {
    because::nimble_harden_samplers(
      mcmc_conf,
      family          = family,
      nimble_samplers = nimble_samplers,
      quiet           = quiet
    )

    if (!quiet) {
      message("Building and compiling NIMBLE MCMC (this may take a moment)...")
    }
    nimble_mcmc <- nimble::buildMCMC(mcmc_conf)
  }
  
  if (!parallel || n.cores == 1 || n.chains == 1) {
    compiled_model <- nimble::compileNimble(nimble_model)
    compiled_mcmc <- nimble::compileNimble(nimble_mcmc, project = nimble_model)
  } else {
    compiled_mcmc <- NULL
    compiled_model <- NULL
  }

  nimble_waic <- NULL

  if (parallel && n.cores > 1 && n.chains > 1) {
    if (!quiet) {
      message(sprintf(
        "Running %d NIMBLE chains in parallel on %d cores...",
        n.chains,
        n.cores
      ))
    }

    if (is.null(cl)) {
      cl <- parallel::makeCluster(n.cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }

    run_nimble_chain <- function(
      chain_id,
      model_string,
      data,
      family,
      nimble_inits,
      monitor,
      n.iter,
      n.burnin,
      n.thin,
      WAIC,
      nimble_samplers,
      quiet
    ) {
      if (!requireNamespace("nimble", quietly = TRUE)) {
        return(NULL)
      }

      nimble_funcs <- list()
      nimble_string <- model_string

      lines <- strsplit(nimble_string, "\n")[[1]]
      lines <- lines[!grepl("logdensity\\.", lines)]
      lines <- lines[!grepl("log_lik_", lines)]
      lines <- lines[!grepl("lik_matrix_", lines)]
      nimble_string <- paste(lines, collapse = "\n")

      for (v in names(family)) {
        fam_obj <- structure(
          list(name = family[[v]]),
          class = c(paste0("because_family_", family[[v]]), "because_family")
        )
        opt_res <- nimble_family_optimization(
          fam_obj,
          nimble_string,
          variable = v
        )
        nimble_string <- opt_res$model_string
        if (length(opt_res$nimble_functions) > 0) {
          nimble_funcs <- c(nimble_funcs, opt_res$nimble_functions)
        }
        if (
          nimble_string != model_string &&
            any(grepl(paste0("z_", v), monitor))
        ) {
          monitor <- setdiff(monitor, paste0("z_", v))
        }
      }
      if (length(nimble_funcs) > 0) {
        unique_names <- unique(names(nimble_funcs))
        for (fn_name in unique_names) {
          assign(fn_name, nimble_funcs[[fn_name]], envir = environment())
          if (startsWith(fn_name, "d")) {
            try(nimble::registerDistributions(nimble_funcs[fn_name]), silent = TRUE)
          }
        }
      }

      nimble_string <- sub("^\\s*model\\s*\\{", "{", nimble_string)
      nimble_code <- parse(text = nimble_string)[[1]]

      curr_inits <- nimble_inits
      if (is.null(curr_inits)) curr_inits <- list()
      
      set.seed(12345 + chain_id)
      for (p_name in names(curr_inits)) {
          val <- curr_inits[[p_name]]
          if (is.numeric(val) && length(val) == 1) {
              if (grepl("^(beta_|alpha_)", p_name)) {
                  curr_inits[[p_name]] <- val + rnorm(1, 0, 0.1)
              } else if (grepl("^(tau_|sigma_)", p_name)) {
                  curr_inits[[p_name]] <- max(0.1, val * exp(rnorm(1, 0, 0.1)))
              }
          }
      }

      requireNamespace("nimble", quietly = TRUE)
      if (!("package:nimble" %in% search())) {
        suppressPackageStartupMessages(attachNamespace("nimble"))
      }
      
      if (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC") {
        requireNamespace("nimbleHMC", quietly = TRUE)
        if (!("package:nimbleHMC" %in% search())) {
          suppressPackageStartupMessages(attachNamespace("nimbleHMC"))
        }
      }

      nimble_constants <- data
      nimble_data <- list()
      if (!is.null(data[["L_multiPhylo"]])) {
          nimble_data[["L_multiPhylo"]] <- data[["L_multiPhylo"]]
          nimble_constants[["L_multiPhylo"]] <- NULL
      }
      if (!is.null(data[["Prec_multiPhylo"]])) {
          nimble_data[["Prec_multiPhylo"]] <- data[["Prec_multiPhylo"]]
          nimble_constants[["Prec_multiPhylo"]] <- NULL
      }
      if (!is.null(data[["L_phylo"]])) {
          nimble_data[["L_phylo"]] <- data[["L_phylo"]]
          nimble_constants[["L_phylo"]] <- NULL
      }

      worker_model <- suppressMessages(suppressWarnings(nimble::nimbleModel(
        code = nimble_code,
        constants = nimble_constants,
        data = nimble_data,
        inits = curr_inits,
        buildDerivs = (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC"),
        calculate = FALSE
      )))

      worker_conf <- nimble::configureMCMC(
        worker_model,
        monitors = monitor,
        enableWAIC = WAIC
      )
      
      if (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC") {
        if (!requireNamespace("nimbleHMC", quietly = TRUE)) {
          stop("The 'nimbleHMC' package is required to use HMC samplers in NIMBLE.")
        }
        worker_mcmc <- nimbleHMC::buildHMC(worker_model)
      } else {
        because::nimble_harden_samplers(
          worker_conf,
          family          = family,
          nimble_samplers = nimble_samplers,
          quiet           = TRUE
        )
        worker_mcmc <- nimble::buildMCMC(worker_conf)
      }
      worker_c_model <- nimble::compileNimble(worker_model)
      worker_c_mcmc <- nimble::compileNimble(
        worker_mcmc,
        project = worker_model
      )

      res <- try({
        samples <- nimble::runMCMC(
          worker_c_mcmc,
          niter = n.iter,
          nburnin = n.burnin,
          nchains = 1,
          thin = n.thin,
          samplesAsCodaMCMC = TRUE,
          WAIC = WAIC
        )
        samples
      }, silent = TRUE)

      if (inherits(res, "try-error")) {
        return(paste("NIMBLE WORKER ERROR:", as.character(res)))
      }
      return(res)
    }

    parallel::clusterExport(
      cl,
      c(
        "model_string", "data", "family", "nimble_inits",
        "monitor", "n.iter", "n.burnin", "n.thin",
        "WAIC", "quiet", "run_nimble_chain", "nimble_samplers"
      ),
      envir = environment()
    )

    if ("package:because" %in% search()) {
      parallel::clusterEvalQ(cl, library(because))
    }
    if ("package:because.phybase" %in% search()) {
      parallel::clusterEvalQ(cl, library(because.phybase))
    }

    chain_results <- parallel::parLapply(cl, seq_len(n.chains), function(i) {
      run_nimble_chain(
        chain_id = i,
        model_string = model_string,
        data = data,
        family = if (!is.null(family)) as.list(family) else NULL,
        nimble_inits = nimble_inits,
        monitor = monitor,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        WAIC = WAIC,
        nimble_samplers = nimble_samplers,
        quiet = quiet
      )
    })

    for (i in seq_along(chain_results)) {
      if (is.character(chain_results[[i]])) {
        stop(sprintf("Chain %d failed: %s", i, chain_results[[i]]))
      }
    }

    if (WAIC && is.list(chain_results[[1]]) && !is.null(chain_results[[1]]$WAIC)) {
        nimble_waic <- chain_results[[1]]$WAIC
    }

    samples <- coda::mcmc.list(lapply(chain_results, function(x) {
        if (is.list(x) && !is.null(x$samples)) x$samples else x
    }))
  } else {
    if (!quiet) {
      message(sprintf(
        "Sampling %d chains sequentially via NIMBLE...",
        n.chains
      ))
    }

    nimble_run_res <- run_nimble_model(
      nimble_code  = nimble_code,
      constants    = data,
      data         = list(),
      inits        = nimble_inits,
      n.chains     = n.chains,
      n.iter       = n.iter,
      n.burnin     = n.burnin,
      n.thin       = n.thin,
      monitor      = monitor,
      quiet        = quiet,
      family       = family,
      nimble_samplers = nimble_samplers,
      nimble_waic  = WAIC
    )
    samples     <- nimble_run_res$samples
    nimble_waic <- nimble_run_res$WAIC
  }
  model <- nimble_model
  saved_nimble_compiled <- if (exists("compiled_mcmc")) compiled_mcmc else NULL
  saved_nimble_cmodel   <- if (exists("compiled_model")) compiled_model else NULL
  saved_nimble_samplers <- nimble_samplers

  return(list(
    samples = samples,
    model = model,
    monitor = monitor,
    nimble_waic = nimble_waic,
    saved_nimble_compiled = saved_nimble_compiled,
    saved_nimble_cmodel = saved_nimble_cmodel,
    saved_nimble_samplers = saved_nimble_samplers
  ))
}
