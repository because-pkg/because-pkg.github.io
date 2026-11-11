#' Build variable-to-level map for Python/NumPyro
#' @keywords internal
.build_py_levels <- function(hierarchical_info) {
  if (is.null(hierarchical_info) || is.null(hierarchical_info$levels)) {
    return(NULL)
  }
  py_levels <- lapply(hierarchical_info$levels, as.character)
  if (!is.null(hierarchical_info$structure_levels)) {
    for (s_name in names(hierarchical_info$structure_levels)) {
      s_lvl <- hierarchical_info$structure_levels[[s_name]]
      if (!is.null(s_lvl) && s_lvl %in% names(py_levels)) {
        py_levels[[s_lvl]] <- unique(c(py_levels[[s_lvl]], s_name))
      }
    }
  }
  if (!is.null(hierarchical_info$link_vars)) {
    for (l_name in names(hierarchical_info$link_vars)) {
      l_var <- hierarchical_info$link_vars[[l_name]]
      if (!is.null(l_var) && l_name %in% names(py_levels)) {
        py_levels[[l_name]] <- unique(c(py_levels[[l_name]], as.character(l_var)))
      }
    }
  }
  lapply(py_levels, as.list)
}

#' Detect Available CPU Cores (Cluster / SLURM / PBS / SGE / LSF / cgroup Safe)
#'
#' @return Integer number of available CPU cores
#' @keywords internal
.detect_available_cores <- function() {
  env_vars <- c(
    "SLURM_CPUS_PER_TASK",
    "SLURM_CPUS_ON_NODE",
    "PBS_NP",
    "NSLOTS",
    "LSB_DJOB_NUMPROC",
    "MC_CORES"
  )
  for (ev in env_vars) {
    val <- trimws(Sys.getenv(ev, unset = ""))
    if (nzchar(val) && grepl("^[0-9]+$", val)) {
      n_ev <- as.integer(val)
      if (!is.na(n_ev) && n_ev > 0L) return(n_ev)
    }
  }
  # Check Linux cgroup v2 / v1 CPU quotas (containers & HPC nodes)
  if (file.exists("/sys/fs/cgroup/cpu.max")) {
    cg2 <- tryCatch(readLines("/sys/fs/cgroup/cpu.max", n = 1L, warn = FALSE), error = function(e) "")
    parts <- strsplit(trimws(cg2), "\\s+")[[1]]
    if (length(parts) == 2L && grepl("^[0-9]+$", parts[1]) && grepl("^[0-9]+$", parts[2])) {
      quota <- as.numeric(parts[1])
      period <- as.numeric(parts[2])
      if (quota > 0 && period > 0) return(max(1L, as.integer(floor(quota / period))))
    }
  }
  if (file.exists("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") && file.exists("/sys/fs/cgroup/cpu/cpu.cfs_period_us")) {
    q_str <- tryCatch(readLines("/sys/fs/cgroup/cpu/cpu.cfs_quota_us", n = 1L, warn = FALSE), error = function(e) "-1")
    p_str <- tryCatch(readLines("/sys/fs/cgroup/cpu/cpu.cfs_period_us", n = 1L, warn = FALSE), error = function(e) "100000")
    quota <- suppressWarnings(as.numeric(trimws(q_str)))
    period <- suppressWarnings(as.numeric(trimws(p_str)))
    if (!is.na(quota) && !is.na(period) && quota > 0 && period > 0) {
      return(max(1L, as.integer(floor(quota / period))))
    }
  }
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.na(detected) || detected < 1L) return(1L)
  max(1L, as.integer(detected - 1L))
}

.ensure_numpyro_ensemble_support <- function(because_py = NULL) {
  tryCatch({
    reticulate::py_run_string("
import os
import sys
import importlib.util
import because.api as _bapi
import because.builder.compiler as _bcomp

_repo_compiler = '/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because_py/because/builder/compiler.py'
if os.path.exists(_repo_compiler) and not hasattr(_bcomp, '_spectral_dcat_enabled'):
    _spec = importlib.util.spec_from_file_location('because.builder.compiler', _repo_compiler)
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    _mod._spectral_dcat_enabled = True
    sys.modules['because.builder.compiler'] = _mod
    _bapi.NumPyroBuilder = _mod.NumPyroBuilder

_repo_api = '/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because_py/because/api.py'
if os.path.exists(_repo_api) and not hasattr(_bapi, '_wave_ensemble_enabled'):
    _was_init = getattr(_bapi, '_jax_initialized', False)
    _spec_a = importlib.util.spec_from_file_location('because.api', _repo_api)
    _mod_a = importlib.util.module_from_spec(_spec_a)
    _mod_a.__package__ = 'because'
    _spec_a.loader.exec_module(_mod_a)
    if _was_init:
        for _attr in ('_jax_initialized', 'jax', 'jnp', 'numpyro', 'dist', 'MCMC', 'NUTS', 'Predictive', 'DiscreteHMCGibbs', 'diag'):
            setattr(_mod_a, _attr, getattr(_bapi, _attr, None))
    _mod_a._wave_ensemble_enabled = True
    sys.modules['because.api'] = _mod_a
    for _k, _v in _mod_a.__dict__.items():
        if not _k.startswith('__'):
            setattr(_bapi, _k, _v)
")
  }, error = function(e) NULL)
}

run_numpyro_model <- function(eq_strings, flat_data, family, priors, py_structures,
                              induced_correlations = NULL,
                              n_chains, n_iter, n_warmup, adapt_delta, max_treedepth,
                              prior_scale_fixed, quiet, WAIC = FALSE, levels = NULL,
                              n_cores = 1L) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required to use engine = 'numpyro'")
  }
  
  because_py <- tryCatch({
    reticulate::import("because.api")
  }, error = function(e) {
    stop("Failed to import 'because.api'. Make sure because_py is installed and reticulate is configured.")
  })
  .ensure_numpyro_ensemble_support(because_py)
  
  py_induced_cors <- if (!is.null(induced_correlations) && length(induced_correlations) > 0) {
    lapply(induced_correlations, as.character)
  } else {
    NULL
  }
  
  # Ensure vector memory limit is not constraining the R process
  tryCatch({
    if (is.finite(mem.maxVSize())) {
      mem.maxVSize(Inf)
    }
  }, error = function(e) NULL)

  py_result <- because_py$fit(
    equations = eq_strings,
    data = flat_data,
    family = family,
    num_chains = as.integer(n_chains),
    num_samples = as.integer(n_iter),
    num_warmup = as.integer(n_warmup),
    n_cores = as.integer(max(1L, n_cores)),
    cor_matrices = py_structures,
    induced_correlations = py_induced_cors,
    adapt_delta = adapt_delta,
    max_treedepth = as.integer(max_treedepth),
    prior_scale_fixed = prior_scale_fixed,
    calculate_waic = WAIC,
    quiet = quiet,
    levels = levels
  )
  
  return(py_result)
}

#' Convert NumPyro output into coda mcmc.list
#'
#' @param py_result Result list from fit_numpyro_model
#' @return coda::mcmc.list
#' @keywords internal
#' @export
format_numpyro_samples <- function(py_result) {
  raw_samples <- py_result$samples
  
  if (length(raw_samples) == 0) return(NULL)
  
  # Determine chains and iterations from the first parameter's dimensions
  first_param <- raw_samples[[1]]
  num_chains <- dim(first_param)[1]
  num_iters <- dim(first_param)[2]
  
  chain_list <- vector("list", num_chains)
  for (ch in 1:num_chains) {
    col_list <- list()
    for (param_name in names(raw_samples)) {
      param_data <- raw_samples[[param_name]]
      p_dims <- dim(param_data)
      p_ndim <- length(p_dims)
      
      if (p_ndim == 2) {
        # Scalar parameter: shape (chains, iters)
        col <- matrix(param_data[ch, ], ncol = 1)
        colnames(col) <- param_name
        col_list[[length(col_list) + 1L]] <- col
      } else if (p_ndim == 3) {
        # Vector parameter: shape (chains, iters, length)
        param_len <- p_dims[3]
        # Only retain small vector parameters like cutpoints (<= 20)
        # Avoid expanding 668-length latent species vectors into hundreds of thousands of columns
        if (param_len <= 20) {
          cols <- matrix(param_data[ch, , ], ncol = param_len)
          colnames(cols) <- paste0(param_name, "[", 1:param_len, "]")
          col_list[[length(col_list) + 1L]] <- cols
        }
      } else if (p_ndim == 4) {
        # Matrix parameter: shape (chains, iters, rows, cols)
        # Specifically extract 2x2 LKJ Cholesky factor L_corr_* -> rho_* if not already present
        if (grepl("^L_corr_", param_name) && p_dims[3] == 2 && p_dims[4] == 2) {
          rho_name <- sub("^L_corr_", "rho_", param_name)
          if (!rho_name %in% names(raw_samples)) {
            col <- matrix(param_data[ch, , 2, 1], ncol = 1)
            colnames(col) <- rho_name
            col_list[[length(col_list) + 1L]] <- col
          }
        }
      }
    }
    chain_mat <- do.call(cbind, col_list)
    chain_list[[ch]] <- coda::mcmc(chain_mat)
  }
  
  return(coda::mcmc.list(chain_list))
}

#' Execute full NumPyro pipeline for because()
#'
#' @keywords internal
run_numpyro_pipeline <- function(
  equations, data, family, structures, priors, random_terms,
  hierarchical_info, is_hierarchical, latent, induced_correlations = NULL, n.chains, n.iter, n.burnin,
  adapt_delta, max_treedepth, prior_scale_fixed, quiet, WAIC,
  model_string, parameter_map, original_call, engine = "numpyro",
  parallel = FALSE, n.cores = 1L
) {
  eq_strings <- sapply(equations, function(eq) paste(deparse(eq), collapse=" "))
  
  # Process structures for NumPyro
  py_structures <- list()
  for (s_name in names(structures)) {
    s_obj <- structures[[s_name]]
    
    # Determine matrix name in the prepared data
    custom_s_name <- get_structure_name_hook(s_obj)
    mat_name <- if (!is.null(custom_s_name)) paste0(custom_s_name, "_", s_name) else paste0("VCV_", s_name)
    is_prec <- FALSE
    if (is.null(data[[mat_name]]) && !is.null(data[[paste0("Prec_", s_name)]])) {
      mat_name <- paste0("Prec_", s_name)
      is_prec <- TRUE
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
          s_tree_method <- attr(s_obj, "tree_method") %||% "dcat"
          s_n_chains <- as.integer(attr(s_obj, "n.chains") %||% 2L)
          s_n_trees <- if (inherits(s_obj, "multiPhylo")) length(s_obj) else 1L
          s_type <- if (inherits(s_obj, "multiPhylo") && identical(s_tree_method, "ensemble")) {
            "multiPhylo_ensemble"
          } else {
            class(s_obj)[1]
          }
          py_structures[[s_name]] <- list(
            matrix = data[[mat_name]],
            transform_func = env[[func_name]],
            type = s_type,
            n_chains = s_n_chains,
            n_trees = as.integer(s_n_trees)
          )
        }
      } else {
        # Fallback to pure covariance matrix
        py_structures[[s_name]] <- if (is_prec) solve(data[[mat_name]]) else data[[mat_name]]
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

  py_levels <- .build_py_levels(hierarchical_info)

  active_cores <- as.integer(if (isTRUE(parallel)) max(1L, n.cores) else 1L)

  py_result <- run_numpyro_model(
    eq_strings = eq_strings,
    flat_data = flat_data,
    family = if (!is.null(family)) as.list(family) else NULL,
    priors = NULL,
    py_structures = py_structures,
    induced_correlations = induced_correlations,
    n_chains = n.chains,
    n_iter = n.iter - n.burnin,
    n_warmup = n.burnin,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    prior_scale_fixed = prior_scale_fixed,
    quiet = quiet,
    WAIC = WAIC,
    levels = py_levels,
    n_cores = active_cores
  )
  mcmc_samples <- format_numpyro_samples(py_result)
  
  # Emergency auto-save checkpoint to disk: guarantees MCMC draws are preserved even if subsequent code fails!
  tryCatch({
    saveRDS(mcmc_samples, file = "because_checkpoint_samples.rds")
  }, error = function(e) NULL)

  result <- list(
    equations = equations,
    family = family,
    latent = latent,
    categorical_vars = attr(data, "categorical_vars"),
    model      = NULL,        # No live model object for NumPyro
    model_code = model_string, # JAGS-equivalent string (kept for reference)
    numpyro_code = if (!is.null(py_result$model_code)) py_result$model_code else NULL,
    parameter_map = parameter_map,
    samples = mcmc_samples,
    data = data,
    original_data = data,
    dsep = NULL,
    priors = priors,
    hierarchical_info = if (!is.null(is_hierarchical) && is_hierarchical) hierarchical_info else NULL,
    induced_correlations = induced_correlations,
    engine = engine,
    parallel = parallel,
    n.cores = n.cores,
    quiet = quiet,
    py_result = list(
      waic = py_result$waic,
      model_code = py_result$model_code,
      model_string = py_result$model_string,
      dsep_results = py_result$dsep_results
    ),
    structures = structures,
    structure = if (!is.null(structures) && length(structures) > 0) structures[[1]] else NULL,
    input = list(
      structure = if (!is.null(structures) && length(structures) > 0) structures[[1]] else NULL,
      structures = structures,
      family = family,
      equations = equations
    )
  )
  if (WAIC && !is.null(py_result$waic)) {
    tryCatch({
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
    }, error = function(e) {
      warning(paste0("Could not construct WAIC table: ", e$message))
    })
  }
  # Compute summary statistics on structural parameters to avoid memory explosion on species-level latent vectors
  pnames_all <- if (!is.null(mcmc_samples)) colnames(mcmc_samples[[1]]) else character(0)
  target_vars <- pnames_all[grepl("^(alpha|beta|sigma|tau|lambda|rho|psi|r_|cutpoint|K_tree)", pnames_all)]
  if (length(target_vars) == 0 && length(pnames_all) > 0) {
    target_vars <- if (length(pnames_all) > 100) head(pnames_all, 100) else pnames_all
  }

  sub_mcmc <- if (length(target_vars) > 0) mcmc_samples[, target_vars, drop = FALSE] else mcmc_samples
  sum_stats <- if (!is.null(sub_mcmc)) summary(sub_mcmc) else NULL

  if (!is.null(sub_mcmc) && length(sub_mcmc) > 1) {
    tryCatch({
      n_ch <- length(sub_mcmc)
      # Convert sub_mcmc chains to numeric matrices ONCE outside loops
      chain_mats <- lapply(sub_mcmc, as.matrix)
      pnames <- colnames(chain_mats[[1]])
      n_params <- length(pnames)
      rhat_vals <- numeric(n_params)
      rhat_max_vals <- numeric(n_params)
      rhat_tree_vals <- numeric(n_params)
      n_iter_chain <- nrow(chain_mats[[1]])

      # Detect multiPhylo ensemble: each chain fixed on a tree and multiple trees represented
      k_col <- grep("^K_tree_", pnames, value = TRUE)[1]
      has_ens_struct <- any(vapply(py_structures, function(ps) is.list(ps) && identical(ps$type, "multiPhylo_ensemble"), logical(1)))
      is_ensemble_run <- FALSE
      chain_trees <- NULL
      trees_multi <- NULL
      if (!is.na(k_col) && has_ens_struct) {
        chain_trees <- vapply(seq_len(n_ch), function(c) {
          as.integer(round(chain_mats[[c]][1, k_col]))
        }, integer(1))
        u_trees <- sort(unique(chain_trees))
        trees_multi <- u_trees[vapply(u_trees, function(tk) sum(chain_trees == tk) >= 2, logical(1))]
        if (length(trees_multi) > 0) {
          is_ensemble_run <- TRUE
        }
      }

      # Precompute chain means and chain variances vectorized across all parameters
      chain_means_mat <- do.call(rbind, lapply(chain_mats, colMeans))
      chain_vars_mat <- do.call(rbind, lapply(chain_mats, function(m) apply(m, 2, stats::var)))

      for (p in 1:n_params) {
        c_means <- chain_means_mat[, p]
        c_vars <- chain_vars_mat[, p]
        B <- n_iter_chain * stats::var(c_means)
        W <- mean(c_vars)
        if (!is.na(W) && W > 1e-12) {
          var_plus <- ((n_iter_chain - 1) / n_iter_chain) * W + (1 / n_iter_chain) * B
          rhat_tree_vals[p] <- sqrt(var_plus / W)
        } else {
          rhat_tree_vals[p] <- 1.0
        }

        if (is_ensemble_run) {
          rw <- vapply(trees_multi, function(tk) {
            idx <- which(chain_trees == tk)
            W_k <- mean(c_vars[idx])
            B_k <- n_iter_chain * stats::var(c_means[idx])
            if (is.na(W_k) || W_k <= 1e-12) return(1.0)
            vp_k <- ((n_iter_chain - 1) / n_iter_chain) * W_k + (1 / n_iter_chain) * B_k
            sqrt(vp_k / W_k)
          }, numeric(1))
          rhat_vals[p] <- mean(rw, na.rm = TRUE)
          rhat_max_vals[p] <- max(rw, na.rm = TRUE)
        } else {
          rhat_vals[p] <- rhat_tree_vals[p]
        }
      }
      if (is_ensemble_run) {
        sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = rhat_vals, Rhat_max = rhat_max_vals)
      } else {
        sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = rhat_vals)
      }
    }, error = function(e) {})
  }
  result$summary <- sum_stats
  mp_structs <- Filter(function(ps) is.list(ps) && ps$type %in% c("multiPhylo", "multiPhylo_ensemble"), py_structures)
  if (length(mp_structs) > 0) {
    result$tree_method <- if (identical(mp_structs[[1]]$type, "multiPhylo_ensemble")) "ensemble" else "dcat"
    result$n_trees_total <- mp_structs[[1]]$n_trees
  }

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

  avail_cores <- .detect_available_cores()
  target_cores <- as.integer(if (parallel) min(n.cores, n.chains) else 1L)
  thread_limit <- as.integer(max(1L, if (parallel) n.cores else avail_cores))
  host_devices <- as.integer(max(target_cores, n.cores, n.chains, avail_cores, 32L))

  .thread_vars <- list(
    OMP_NUM_THREADS            = "1",
    OPENBLAS_NUM_THREADS       = "1",
    GOTO_NUM_THREADS           = "1",
    MKL_NUM_THREADS            = "1",
    MKL_DOMAIN_NUM_THREADS     = "1",
    NUMEXPR_NUM_THREADS        = "1",
    LLVM_NUM_THREADS           = "1",
    TF_NUM_INTEROP_THREADS     = as.character(thread_limit),
    TF_NUM_INTRAOP_THREADS     = as.character(thread_limit),
    XLA_PYTHON_CLIENT_PREALLOCATE = "false"
  )

  for (.v in names(.thread_vars)) {
    if (!nzchar(Sys.getenv(.v))) {
      do.call(Sys.setenv, stats::setNames(list(.thread_vars[[.v]]), .v))
    }
  }

  .current_xla <- Sys.getenv("XLA_FLAGS")
  if (grepl("--xla_force_host_platform_device_count=[0-9]+", .current_xla)) {
    .current_xla <- sub(
      "--xla_force_host_platform_device_count=[0-9]+",
      paste0("--xla_force_host_platform_device_count=", host_devices),
      .current_xla
    )
  } else {
    .current_xla <- trimws(paste(.current_xla, paste0("--xla_force_host_platform_device_count=", host_devices)))
  }
  .xla_additions <- character(0)
  if (!grepl("--xla_cpu_multi_thread_eigen", .current_xla))
    .xla_additions <- c(.xla_additions, "--xla_cpu_multi_thread_eigen=false")
  if (!grepl("intra_op_parallelism_threads", .current_xla))
    .xla_additions <- c(.xla_additions, paste0("intra_op_parallelism_threads=", thread_limit))
  if (!grepl("inter_op_parallelism_threads", .current_xla))
    .xla_additions <- c(.xla_additions, paste0("inter_op_parallelism_threads=", thread_limit))
  Sys.setenv(XLA_FLAGS = trimws(paste(.current_xla, paste(.xla_additions, collapse = " "))))

  b_api <- tryCatch({
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

  # Ensure XLA_FLAGS and _init_jax(host_devices) run BEFORE any structure py_code imports jax
  .py_set_code <- paste(
    "import os, sys",
    paste(sapply(names(.thread_vars), function(.v) {
      sprintf("os.environ.setdefault('%s', '%s')", .v, .thread_vars[[.v]])
    }), collapse = "\n"),
    sprintf("if 'jax' not in sys.modules:\n    os.environ['XLA_FLAGS'] = '%s'", Sys.getenv("XLA_FLAGS")),
    sep = "\n"
  )
  tryCatch({
    reticulate::py_run_string(.py_set_code)
    b_api$`_init_jax`(as.integer(host_devices))
  }, error = function(e) NULL)

  return(b_api)
}

