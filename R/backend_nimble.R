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
