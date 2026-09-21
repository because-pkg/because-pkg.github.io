#' Run a NumPyro model via because_py
#'
#' @param eq_strings Character vector of equations
#' @param flat_data Flattened list of data
#' @param family Named character vector of families
#' @param priors Named list of priors
#' @param py_structures Processed structural matrices/functions for NumPyro
#' @param n_chains Number of MCMC chains
#' @param n_iter Number of iterations
#' @param n_warmup Number of warmup iterations
#' @param adapt_delta Target acceptance probability
#' @param max_treedepth Maximum tree depth
#' @param prior_scale_fixed Scale factor for fixed effects
#' @param quiet Logical, suppress output
#' @return Raw result from because_py$fit_numpyro_model
#' @export
run_numpyro_model <- function(eq_strings, flat_data, family, priors, py_structures,
                              n_chains, n_iter, n_warmup, adapt_delta, max_treedepth,
                              prior_scale_fixed, quiet) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required to use engine = 'numpyro'")
  }
  
  because_py <- tryCatch({
    reticulate::import("because.api")
  }, error = function(e) {
    stop("Failed to import 'because.api'. Make sure because_py is installed and reticulate is configured.")
  })
  
  py_result <- because_py$fit_numpyro_model(
    equations = eq_strings,
    data = flat_data,
    family = family,
    priors = priors,
    n_chains = as.integer(n_chains),
    n_samples = as.integer(n_iter),
    n_warmup = as.integer(n_warmup),
    cor_matrices = py_structures,
    adapt_delta = adapt_delta,
    max_treedepth = as.integer(max_treedepth),
    prior_scale_fixed = prior_scale_fixed
  )
  
  return(py_result)
}

#' Convert NumPyro output into coda mcmc.list
#'
#' @param py_result Result list from fit_numpyro_model
#' @return coda::mcmc.list
#' @export
format_numpyro_samples <- function(py_result) {
  raw_samples <- py_result$samples
  
  if (length(raw_samples) == 0) return(NULL)
  
  # Determine chains and iterations from the first parameter's dimensions
  first_param <- raw_samples[[1]]
  num_chains <- dim(first_param)[1]
  num_iters <- dim(first_param)[2]
  
  chain_list <- list()
  for (ch in 1:num_chains) {
    chain_mat <- NULL
    for (param_name in names(raw_samples)) {
      param_data <- raw_samples[[param_name]]
      
      if (length(dim(param_data)) == 2) {
        # Scalar parameter: shape (chains, iters)
        col <- matrix(param_data[ch, ], ncol = 1)
        colnames(col) <- param_name
        chain_mat <- if (is.null(chain_mat)) col else cbind(chain_mat, col)
      } else if (length(dim(param_data)) == 3) {
        # Vector parameter: shape (chains, iters, length)
        param_len <- dim(param_data)[3]
        cols <- matrix(param_data[ch, , ], ncol = param_len)
        colnames(cols) <- paste0(param_name, "[", 1:param_len, "]")
        chain_mat <- if (is.null(chain_mat)) cols else cbind(chain_mat, cols)
      }
    }
    chain_list[[ch]] <- coda::mcmc(chain_mat)
  }
  
  return(coda::mcmc.list(chain_list))
}

#' Execute full NumPyro pipeline for because()
#'
#' @keywords internal
run_numpyro_pipeline <- function(
  equations, data, family, structures, priors, random_terms,
  hierarchical_info, is_hierarchical, latent, n.chains, n.iter, n.burnin,
  adapt_delta, max_treedepth, prior_scale_fixed, quiet, WAIC,
  model_string, parameter_map, original_call, engine = "numpyro"
) {
  eq_strings <- sapply(equations, function(eq) paste(deparse(eq), collapse=" "))
  
  # Process structures for NumPyro
  py_structures <- list()
  for (s_name in names(structures)) {
    s_obj <- structures[[s_name]]
    
    # Determine matrix name in the prepared data
    custom_s_name <- get_structure_name_hook(s_obj)
    mat_name <- if (!is.null(custom_s_name)) paste0(custom_s_name, "_", s_name) else paste0("VCV_", s_name)
    if (is.null(data[[mat_name]]) && !is.null(data[[paste0("Prec_", s_name)]])) {
      mat_name <- paste0("Prec_", s_name)
    }
    
    if (!is.null(data[[mat_name]])) {
      # Ask the extension package for the Python JAX code
      py_code <- numpyro_structure_definition(s_obj, engine = "numpyro")
      if (!is.null(py_code)) {
        # Compile into a Python function using reticulate
        env <- reticulate::py_run_string(py_code)
        funcs <- names(env)
        expected_name <- paste0(s_name, "_transform")
        func_name <- if (expected_name %in% funcs) expected_name else funcs[!funcs %in% c("numpyro", "jnp", "jax", "dist", "np", "r")][1]
        if (!is.null(func_name) && (func_name %in% names(env))) {
          py_structures[[s_name]] <- list(
            matrix = data[[mat_name]],
            transform_func = env[[func_name]],
            type = class(s_obj)[1]
          )
        }
      } else {
        # Fallback to pure matrix
        py_structures[[s_name]] <- data[[mat_name]]
      }
      
      # Append + (1 | s_name) to valid endogenous equations
      for (i in seq_along(equations)) {
         eq <- equations[[i]]
         response <- trimws(strsplit(deparse(eq), "~")[[1]][1])
         
         # Apply if dimension matches or valid level
         is_valid <- TRUE
         if (!is.null(is_hierarchical) && is_hierarchical && !is.null(hierarchical_info)) {
             s_lvl <- hierarchical_info$structure_levels[[s_name]]
             tryCatch({
                 resp_lvl <- infer_variable_level(response, hierarchical_info$levels, data = NULL, equations = equations, latent = latent, hierarchy = hierarchical_info$hierarchy)
                 if (is.null(resp_lvl)) {
                     is_valid <- FALSE
                 } else if (!is_valid_structure_mapping_dsep(s_lvl, resp_lvl, hierarchical_info)) {
                     is_valid <- FALSE
                 }
             }, error = function(e) {
                 is_valid <<- FALSE
             })
         }
         
         if (is_valid) {
             eq_str <- eq_strings[i]
             if (!grepl(paste0("\\(1\\s*\\|\\s*", s_name, "\\)"), eq_str)) {
                 eq_strings[i] <- paste0(eq_str, " + (1|", s_name, ")")
             }
         }
      }
    }
  }
  
  flat_data <- flatten_for_python(data)
  
  # Ensure zero-indexing for Python categorical variables
  for (s_name in names(structures)) {
      if (is.null(flat_data[[s_name]])) {
          N_val <- NULL
          if (!is.null(is_hierarchical) && is_hierarchical && !is.null(hierarchical_info)) {
              s_lvl <- hierarchical_info$structure_levels[[s_name]]
              if (!is.null(s_lvl)) {
                  N_val <- flat_data[[paste0("N_", s_lvl)]]
              }
          }
          if (is.null(N_val)) {
              N_val <- if (!is.null(flat_data[["N"]])) flat_data[["N"]] else length(flat_data[[1]])
          }
          flat_data[[s_name]] <- 0:(N_val - 1)
      }
  }
  idx_vars <- grep("_idx", names(flat_data), value = TRUE)
  idx_vars <- c(idx_vars, names(structures))
  for (eq in equations) {
      eq_str <- if (is.character(eq)) eq else paste(deparse(eq), collapse=" ")
      matches <- regmatches(eq_str, gregexpr("\\(1\\s*\\|\\s*[^)]+\\)", eq_str))[[1]]
      for (m in matches) {
          grp <- trimws(strsplit(m, "\\|")[[1]][2])
          grp <- gsub("\\)", "", grp)
          idx_vars <- c(idx_vars, grp)
      }
  }
  idx_vars <- unique(idx_vars)
  for (s_name in idx_vars) {
      if (s_name %in% names(flat_data) && min(flat_data[[s_name]], na.rm=TRUE) >= 1) {
          flat_data[[s_name]] <- as.integer(flat_data[[s_name]] - 1L)
      }
  }
  
  # Re-inject standard random terms that were stripped for JAGS
  if (length(random_terms) > 0) {
      for (rt in random_terms) {
          for (i in seq_along(equations)) {
              resp <- trimws(strsplit(deparse(equations[[i]]), "~")[[1]][1])
              if (resp == rt$response) {
                  re_str <- paste0("\\(1\\s*\\|\\s*", rt$group, "\\)")
                  if (!grepl(re_str, eq_strings[i])) {
                      eq_strings[i] <- paste0(eq_strings[i], " + (1|", rt$group, ")")
                  }
              }
          }
      }
  }

  py_result <- run_numpyro_model(
    eq_strings = eq_strings,
    flat_data = flat_data,
    family = if (!is.null(family)) as.list(family) else NULL,
    priors = NULL,
    py_structures = py_structures,
    n_chains = n.chains,
    n_iter = n.iter - n.burnin,
    n_warmup = n.burnin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    prior_scale_fixed = prior_scale_fixed,
    quiet = quiet
  )
  mcmc_samples <- format_numpyro_samples(py_result)
  
  result <- list(
    equations = equations,
    model      = NULL,        # No live model object for NumPyro
    model_code = model_string, # JAGS-equivalent string (kept for reference)
    numpyro_code = if (!is.null(py_result$model_code)) py_result$model_code else NULL,
    parameter_map = parameter_map,
    samples = mcmc_samples,
    data = data,
    dsep = NULL,
    priors = priors,
    hierarchical_info = if (!is.null(is_hierarchical) && is_hierarchical) hierarchical_info else NULL,
    engine = engine,
    quiet = quiet
  )
  if (WAIC && !is.null(py_result$waic)) {
    waic_df <- data.frame(
      Estimate = c(py_result$waic$elpd_waic$Estimate, py_result$waic$p_waic$Estimate, py_result$waic$waic$Estimate),
      SE = c(py_result$waic$elpd_waic$SE, py_result$waic$p_waic$SE, py_result$waic$waic$SE),
      row.names = c("elpd_waic", "p_waic", "waic")
    )
    
    attr(waic_df, "pointwise") <- data.frame(
      elpd_waic = as.numeric(py_result$waic$pointwise$elpd_waic_i),
      p_waic = as.numeric(py_result$waic$pointwise$p_waic_i),
      waic = as.numeric(py_result$waic$pointwise$waic_i)
    )
    
    attr(waic_df, "dims") <- c(n_obs = py_result$waic$n_obs, n_samples = py_result$waic$n_samples)
    class(waic_df) <- c("because_waic", "data.frame")
    
    result$WAIC <- waic_df
  }
  # Compute summary statistics
  sum_stats <- if (!is.null(mcmc_samples)) summary(mcmc_samples) else NULL
  if (!is.null(mcmc_samples) && n.chains > 1) {
    tryCatch({
      n_ch <- length(mcmc_samples)
      first_chain <- as.matrix(mcmc_samples[[1]])
      pnames <- colnames(first_chain)
      n_params <- length(pnames)
      rhat_vals <- numeric(n_params)
      n_iter_chain <- nrow(first_chain)
      for (p in 1:n_params) {
        chain_means <- numeric(n_ch)
        chain_vars <- numeric(n_ch)
        for (c in 1:n_ch) {
          vals <- as.matrix(mcmc_samples[[c]])[, p]
          chain_means[c] <- mean(vals)
          chain_vars[c] <- var(vals)
        }
        grand_mean <- mean(chain_means)
        B <- n_iter_chain * var(chain_means)
        W <- mean(chain_vars)
        if (W > 0) {
          var_plus <- ((n_iter_chain - 1) / n_iter_chain) * W + (1 / n_iter_chain) * B
          rhat_vals[p] <- sqrt(var_plus / W)
        } else {
          rhat_vals[p] <- 1.0
        }
      }
      sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = rhat_vals)
    }, error = function(e) {})
  }
  result$summary <- sum_stats

  result$call <- original_call
  class(result) <- "because"
  return(result)
}


#' Setup Python/NumPyro Environment and Thread Constraints
#'
#' @param parallel Logical, whether chains run in parallel
#' @param n.cores Number of CPU cores requested
#' @param n.chains Number of MCMC chains
#' @return imported because.api module
#' @noRd
setup_numpyro_environment <- function(parallel = FALSE, n.cores = 1, n.chains = 3) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required when engine = 'numpyro'.")
  }

  tryCatch({
    if (reticulate::virtualenv_exists("because_env")) {
      reticulate::use_virtualenv("because_env", required = TRUE)
    } else if (reticulate::condaenv_exists("because_env")) {
      reticulate::use_condaenv("because_env", required = TRUE)
    }
  }, error = function(e) NULL)

  target_cores <- as.integer(if (parallel) min(n.cores, n.chains) else 1L)

  .thread_vars <- list(
    OMP_NUM_THREADS            = "1",
    OPENBLAS_NUM_THREADS       = "1",
    GOTO_NUM_THREADS           = "1",
    MKL_NUM_THREADS            = "1",
    MKL_DOMAIN_NUM_THREADS     = "1",
    NUMEXPR_NUM_THREADS        = "1",
    LLVM_NUM_THREADS           = "1",
    TF_NUM_INTEROP_THREADS     = as.character(target_cores),
    TF_NUM_INTRAOP_THREADS     = as.character(target_cores),
    XLA_PYTHON_CLIENT_PREALLOCATE = "false"
  )

  for (.v in names(.thread_vars)) {
    if (!nzchar(Sys.getenv(.v))) {
      do.call(Sys.setenv, stats::setNames(list(.thread_vars[[.v]]), .v))
    }
  }

  .current_xla <- Sys.getenv("XLA_FLAGS")
  .xla_additions <- character(0)
  if (!grepl("--xla_force_host_platform_device_count", .current_xla))
    .xla_additions <- c(.xla_additions,
                        paste0("--xla_force_host_platform_device_count=", target_cores))
  if (!grepl("--xla_cpu_multi_thread_eigen", .current_xla))
    .xla_additions <- c(.xla_additions, "--xla_cpu_multi_thread_eigen=false")
  if (!grepl("intra_op_parallelism_threads", .current_xla))
    .xla_additions <- c(.xla_additions, paste0("intra_op_parallelism_threads=", target_cores))
  if (!grepl("inter_op_parallelism_threads", .current_xla))
    .xla_additions <- c(.xla_additions, paste0("inter_op_parallelism_threads=", target_cores))
  if (length(.xla_additions) > 0)
    Sys.setenv(XLA_FLAGS = trimws(paste(.current_xla, paste(.xla_additions, collapse = " "))))

  if (reticulate::py_available(initialize = FALSE)) {
    .py_set_code <- paste(
      "import os",
      paste(sapply(names(.thread_vars), function(.v) {
        sprintf("os.environ.setdefault('%s', '%s')", .v, .thread_vars[[.v]])
      }), collapse = "\n"),
      sprintf("os.environ.setdefault('XLA_FLAGS', '%s')", Sys.getenv("XLA_FLAGS")),
      sep = "\n"
    )
    tryCatch(
      reticulate::py_run_string(.py_set_code),
      error = function(e) NULL
    )
  }

  tryCatch({
    reticulate::import("because.api")
  }, error = function(e) {
    current_env <- "unknown"
    tryCatch({
      current_env <- reticulate::py_config()$python
    }, error = function(e) {})
    stop(sprintf(paste0(
      "Failed to import python module 'because.api'.\n",
      "Python is currently running from: %s\n\n",
      "This usually means Python was initialized to a different environment\n",
      "before 'library(because)' was called (e.g. by RStudio or another package).\n\n",
      "Quick fix --- add this line BEFORE library(because) in your script:\n",
      "  reticulate::use_virtualenv('because_env', required = TRUE)\n\n",
      "If because_env does not exist yet, install it first with:\n",
      "  install_because_numpyro()"
    ), current_env))
  })
}

