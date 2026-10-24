#' Run d-separation tests using NumPyro backend
#'
#' Evaluates the basis set of conditional independence tests using NumPyro,
#' handling structure transforms, Python data flattening, cross-scale tests,
#' and caching.
#'
#' @keywords internal
run_dsep_numpyro <- function(
  equations,
  structures,
  data,
  latent,
  random_terms,
  hierarchical_info,
  all_poly_terms,
  family,
  quiet,
  original_call,
  original_data,
  id_col = NULL,
  variability = NULL,
  adapt_delta = 0.95,
  max_treedepth = 10,
  expand_ordered = FALSE,
  reuse_models = NULL,
  aggregate_crossscale = NULL,
  n.iter = 2000,
  n.burnin = 1000,
  n.chains = 4,
  n.thin = 1,
  parallel = FALSE,
  n.cores = 1,
  dsep_max_obs = 10000,
  prior_scale_fixed = NULL,
  engine = "numpyro",
  structure = NULL,
  because_py = NULL
) {
  if (is.null(because_py)) {
    because_py <- tryCatch({
      reticulate::import("because.api", delay_load = FALSE)
    }, error = function(e) {
      stop("Failed to import 'because.api'. Make sure because_py is installed and reticulate is configured.")
    })
  }

  eq_strings <- sapply(equations, function(eq) paste(deparse(eq), collapse = " "))
  
  # Process structures for NumPyro (early evaluation for dsep)
  py_structures <- list()
  for (s_name in names(structures)) {
    s_obj <- structures[[s_name]]
    custom_s_name <- get_structure_name_hook(s_obj)
    mat_name <- if (!is.null(custom_s_name)) paste0(custom_s_name, "_", s_name) else paste0("VCV_", s_name)
    if (is.null(data[[mat_name]]) && !is.null(data[[paste0("Prec_", s_name)]])) {
      mat_name <- paste0("Prec_", s_name)
    }
    if (!is.null(data[[mat_name]])) {
      py_code <- numpyro_structure_definition(s_obj, engine = "numpyro")
      if (!is.null(py_code)) {
        env <- reticulate::py_run_string(py_code)
        funcs <- names(env)
        expected_name <- paste0(s_name, "_transform")
        func_name <- if (expected_name %in% funcs) expected_name else funcs[!funcs %in% c("numpyro", "jnp", "jax", "dist", "np", "r")][1]
        if (!is.null(func_name) && (func_name %in% names(env))) {
          py_structures[[s_name]] <- list(matrix = data[[mat_name]], transform_func = env[[func_name]], type = class(s_obj)[1])
        }
      } else {
        py_structures[[s_name]] <- data[[mat_name]]
      }
    }
  }
  
  # --- HOTFIX: Ultimate safety net for Python/JAX ---
  # Remove any character or factor vectors from data before passing to NumPyro
  # to guarantee we never hit a JAX string dtype error.
  for (n in names(data)) {
    if (is.character(data[[n]]) || is.factor(data[[n]])) {
      data[[n]] <- NULL
    }
  }
  
  flat_data <- flatten_for_python(data)
  for (s_name in names(structures)) {
    if (is.null(flat_data[[s_name]])) {
      N_val <- if (!is.null(flat_data[["N"]])) flat_data[["N"]] else length(flat_data[[1]])
      flat_data[[s_name]] <- 0:(N_val - 1)
    }
  }
  idx_vars <- grep("_idx", names(flat_data), value = TRUE)
  idx_vars <- c(idx_vars, names(structures))
  for (eq in equations) {
    eq_str <- if (is.character(eq)) eq else paste(deparse(eq), collapse = " ")
    matches <- regmatches(eq_str, gregexpr("\\(1\\s*\\|\\s*[^)]+\\)", eq_str))[[1]]
    for (m in matches) {
      grp <- trimws(strsplit(m, "\\|")[[1]][2])
      grp <- gsub("\\)", "", grp)
      idx_vars <- c(idx_vars, grp)
    }
  }
  idx_vars <- unique(idx_vars)
  for (s_name in idx_vars) {
    if (s_name %in% names(flat_data) && min(flat_data[[s_name]], na.rm = TRUE) >= 1) {
      flat_data[[s_name]] <- as.integer(flat_data[[s_name]] - 1L)
    }
  }

  # Compute the optimal DAG/MAG tests natively using because_dsep (with dagitty)
  native_dsep <- because_dsep(
    equations = equations,
    latent = latent,
    random_terms = random_terms,
    hierarchical_info = hierarchical_info,
    poly_terms = all_poly_terms,
    categorical_vars = attr(data, "categorical_vars"),
    family = family,
    quiet = TRUE
  )
  
  current_tests <- list()
  if (!is.null(native_dsep)) {
    if (!is.null(latent) && length(latent) > 0) {
      current_tests <- native_dsep$tests
    } else {
      current_tests <- native_dsep
    }
  }
  
  # --- A priori cross-hierarchy filter ---
  # Tests where the response and the focal predictor belong to orthogonal
  # hierarchical branches are trivially satisfied. Skip them.
  if (length(current_tests) > 0) {
    cross_hier_flags <- sapply(current_tests, function(eq) {
      is_cross_hierarchy_test(eq, hierarchical_info)
    })
    n_cross <- sum(cross_hier_flags)
    if (n_cross > 0 && !quiet) {
      skipped_labels <- sapply(current_tests[cross_hier_flags], function(eq) {
        resp     <- as.character(eq)[2]
        test_var <- attr(eq, "test_var")
        if (!is.null(test_var)) {
          sprintf("  %s _||_ %s (orthogonal hierarchy branches - trivially satisfied)",
                  resp, test_var)
        } else {
          paste(deparse(eq), collapse = " ")
        }
      })
      message(sprintf(
        "\nSkipping %d cross-hierarchy d-sep test(s) (trivially satisfied by design):",
        n_cross
      ))
      for (lbl in skipped_labels) message(lbl)
    }
    current_tests <- current_tests[!cross_hier_flags]
  }
  
  incremental_check <- find_reusable_tests(
    current_tests,
    equations,
    reuse_models,
    data,
    family = if (!is.null(family)) as.list(family) else NULL,
    quiet = quiet
  )
  
  reused_results <- incremental_check$found
  tests_to_run_indices <- incremental_check$missing_indices
  
  dsep_equations_to_run <- NULL
  if (length(tests_to_run_indices) == 0 && length(current_tests) > 0) {
    # All tests cached!
    py_result <- list(dsep_results = list())
  } else {
    if (length(tests_to_run_indices) > 0) {
      # Format tests into structured list for Python to bypass its internal graph logic
      dsep_equations_to_run <- list()
      cs_results <- list()
      
      for (i in tests_to_run_indices) {
        eq <- current_tests[[i]]
        
        cs_info <- detect_crossscale_dsep(eq, hierarchical_info)
        do_aggregate <- FALSE
        if (cs_info$is_crossscale) {
          if (is.character(aggregate_crossscale) && length(aggregate_crossscale) == 1 && aggregate_crossscale == "all") {
            do_aggregate <- TRUE
          } else if (is.numeric(aggregate_crossscale) && i %in% aggregate_crossscale) {
            do_aggregate <- TRUE
          }
        }

        if (do_aggregate) {
          if (!quiet) {
            message(sprintf("  -> Cross-scale test detected: focal predictor '%s' at %s level, response '%s' at %s level", cs_info$test_var, cs_info$predictor_level, cs_info$response, cs_info$response_level))
          }
          # Run PGLS in R instead of passing to Python
          synth_iter <- max(1000, n.iter)
          mcmc_res <- run_crossscale_dsep_pgls(
            i = i, test_eq = eq, cs_info = cs_info,
            original_data = original_data, hierarchical_info = hierarchical_info,
            structure = structure, family = family, engine = engine, 
            n.iter = synth_iter, n.chains = n.chains, quiet = quiet
          )
          
          samples_mcmc <- mcmc_res$samples
          param_name <- colnames(as.matrix(samples_mcmc))[1]
          vec <- as.numeric(as.matrix(samples_mcmc)[, 1])
          
          r_hat_val <- 1.0
          n_eff_val <- length(vec)
          if (inherits(samples_mcmc, "mcmc.list") && length(samples_mcmc) > 1) {
            tryCatch({
              n_ch <- length(samples_mcmc)
              chain_means <- numeric(n_ch)
              chain_vars <- numeric(n_ch)
              n_iter_ch <- nrow(as.matrix(samples_mcmc[[1]]))
              for (c in 1:n_ch) {
                vals <- as.numeric(as.matrix(samples_mcmc[[c]])[, 1])
                chain_means[c] <- mean(vals, na.rm = TRUE)
                chain_vars[c] <- var(vals, na.rm = TRUE)
              }
              B <- n_iter_ch * var(chain_means, na.rm = TRUE)
              W <- mean(chain_vars, na.rm = TRUE)
              if (W > 0) {
                var_plus <- ((n_iter_ch - 1) / n_iter_ch) * W + (1 / n_iter_ch) * B
                r_hat_val <- sqrt(var_plus / W)
              }
              
              if (requireNamespace("coda", quietly = TRUE)) {
                n_eff_val <- coda::effectiveSize(samples_mcmc[, param_name, drop = FALSE])
              }
            }, error = function(e) {})
          }
          
          eq_full <- current_tests[[i]]
          resp_name <- as.character(eq_full)[2]
          test_var <- attr(eq_full, "test_var")
          if (is.null(test_var)) {
            rhs <- labels(stats::terms(eq_full))
            if (length(rhs) > 0) test_var <- rhs[1]
          }
          formula_str <- paste(deparse(eq_full), collapse = " ")
          rhs_full <- sub("^[^~]+~\\s*", "", formula_str)
          all_terms <- trimws(strsplit(rhs_full, "\\+")[[1]])
          cond_terms <- all_terms[all_terms != test_var]

          test_str_base <- paste0(
            resp_name,
            " _||_ ",
            test_var,
            if (length(cond_terms) > 0) {
              paste0(" | {", paste(cond_terms, collapse = ","), "}")
            } else {
              " | {} "
            }
          )
          
          test_scale <- attr(eq_full, "scale")
          if (!is.null(test_scale)) {
            test_str_base <- paste0(test_str_base, " [Scale: ", test_scale, "]")
          }
          
          cs_df <- data.frame(
            claim = test_str_base,
            coefficient = param_name,
            mean = round(mean(vec, na.rm = TRUE), 3),
            ci_2.5 = round(unname(quantile(vec, 0.025, na.rm = TRUE)), 3),
            ci_97.5 = round(unname(quantile(vec, 0.975, na.rm = TRUE)), 3),
            rhat = round(r_hat_val, 3),
            n_eff = round(n_eff_val, 0),
            stringsAsFactors = FALSE
          )
          cs_results[[as.character(i)]] <- cs_df
        } else {
          # Standard test, prepare for Python
          resp <- as.character(eq)[2]
          test_var <- attr(eq, "test_var")
          
          rhs_str <- paste(deparse(eq[[3]]), collapse = " ")
          rhs_parts <- trimws(strsplit(rhs_str, "\\+")[[1]])
          cond_set <- rhs_parts[rhs_parts != test_var]
          
          if (is.null(test_var)) {
            vars <- all.vars(eq)
            test_var <- setdiff(vars, resp)[1]
            cond_set <- setdiff(vars, c(resp, test_var))
          }
          
          dsep_equations_to_run[[length(dsep_equations_to_run) + 1]] <- list(
            type = "dsep",
            response = resp,
            test_node = test_var,
            conditioning_set = as.list(cond_set),
            equation_string = paste(deparse(eq), collapse = " ")
          )
        }
      }
    }
    
    py_result <- because_py$fit(
      equations = eq_strings,
      data = flat_data,
      family = if (!is.null(family)) as.list(family) else NULL,
      latent = latent,
      dsep = TRUE,
      dsep_only = TRUE,
      calculate_waic = FALSE,
      num_samples = as.integer(n.iter - n.burnin),
      num_warmup = as.integer(n.burnin),
      num_chains = as.integer(n.chains),
      thinning = as.integer(n.thin),
      n_cores = as.integer(if (parallel) min(n.cores, n.chains) else 1),
      dsep_max_obs = as.integer(dsep_max_obs),
      quiet = quiet,
      cor_matrices = py_structures,
      dsep_equations_to_run = if (!is.null(dsep_equations_to_run)) as.list(dsep_equations_to_run) else NULL,
      prior_scale_fixed = prior_scale_fixed
    )
  }

  new_dsep_df <- NULL
  if (!is.null(py_result$dsep_results) && length(py_result$dsep_results) > 0) {
    new_dsep_df <- do.call(rbind, lapply(py_result$dsep_results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
    if (nrow(new_dsep_df) > 0) {
      names(new_dsep_df)[names(new_dsep_df) == "claim"] <- "Test"
      names(new_dsep_df)[names(new_dsep_df) == "coefficient"] <- "Parameter"
      names(new_dsep_df)[names(new_dsep_df) == "mean"] <- "Estimate"
      names(new_dsep_df)[names(new_dsep_df) == "ci_2.5"] <- "LowerCI"
      names(new_dsep_df)[names(new_dsep_df) == "ci_97.5"] <- "UpperCI"
      names(new_dsep_df)[names(new_dsep_df) == "rhat"] <- "Rhat"
      names(new_dsep_df)[names(new_dsep_df) == "n_eff"] <- "n.eff"
      new_dsep_df$Scale <- NA
      new_dsep_df$is_independent <- NULL
      new_dsep_df$equation <- NULL
    }
  }
  
  # Merge reused results back into the final dataframe in the original order
  dsep_df_list <- list()
  new_idx <- 1
  if (length(current_tests) > 0) {
    for (i in seq_along(current_tests)) {
      if (!is.null(reused_results) && length(reused_results) >= i && !is.null(reused_results[[i]])) {
        dsep_df_list[[i]] <- reused_results[[i]]
      } else if (!is.null(cs_results[[as.character(i)]])) {
        dsep_df_list[[i]] <- cs_results[[as.character(i)]]
        names(dsep_df_list[[i]])[names(dsep_df_list[[i]]) == "claim"] <- "Test"
        names(dsep_df_list[[i]])[names(dsep_df_list[[i]]) == "coefficient"] <- "Parameter"
        names(dsep_df_list[[i]])[names(dsep_df_list[[i]]) == "mean"] <- "Estimate"
        names(dsep_df_list[[i]])[names(dsep_df_list[[i]]) == "ci_2.5"] <- "LowerCI"
        names(dsep_df_list[[i]])[names(dsep_df_list[[i]]) == "ci_97.5"] <- "UpperCI"
        names(dsep_df_list[[i]])[names(dsep_df_list[[i]]) == "rhat"] <- "Rhat"
        names(dsep_df_list[[i]])[names(dsep_df_list[[i]]) == "n_eff"] <- "n.eff"
        dsep_df_list[[i]]$Scale <- NA
        dsep_df_list[[i]]$is_independent <- NULL
        dsep_df_list[[i]]$equation <- NULL
      } else if (!is.null(new_dsep_df) && new_idx <= nrow(new_dsep_df)) {
        dsep_df_list[[i]] <- new_dsep_df[new_idx, , drop = FALSE]
        new_idx <- new_idx + 1
      }

      # Force uniform JAGS-style formatting for all NumPyro results 
      # (including reused results cached from older versions)
      if (!is.null(dsep_df_list[[i]])) {
        eq_full <- current_tests[[i]]
        resp_name <- as.character(eq_full)[2]
        test_var <- attr(eq_full, "test_var")
        if (is.null(test_var)) {
          rhs <- labels(stats::terms(eq_full))
          if (length(rhs) > 0) test_var <- rhs[1]
        }
        formula_str <- paste(deparse(eq_full), collapse = " ")
        rhs_full <- sub("^[^~]+~\\s*", "", formula_str)
        all_terms <- trimws(strsplit(rhs_full, "\\+")[[1]])
        cond_terms <- all_terms[all_terms != test_var]

        test_str_base <- paste0(
          resp_name,
          " _||_ ",
          test_var,
          if (length(cond_terms) > 0) {
            paste0(" | {", paste(cond_terms, collapse = ","), "}")
          } else {
            " | {} "
          }
        )
        
        test_scale <- attr(eq_full, "scale")
        if (!is.null(test_scale)) {
          test_str_base <- paste0(test_str_base, " [Scale: ", test_scale, "]")
        }
        
        dsep_df_list[[i]]$Test <- test_str_base

        # Also force rounding for standard Python tests so they look beautiful
        if ("Estimate" %in% names(dsep_df_list[[i]])) dsep_df_list[[i]]$Estimate <- round(as.numeric(dsep_df_list[[i]]$Estimate), 3)
        if ("LowerCI" %in% names(dsep_df_list[[i]])) dsep_df_list[[i]]$LowerCI <- round(as.numeric(dsep_df_list[[i]]$LowerCI), 3)
        if ("UpperCI" %in% names(dsep_df_list[[i]])) dsep_df_list[[i]]$UpperCI <- round(as.numeric(dsep_df_list[[i]]$UpperCI), 3)
        if ("Rhat" %in% names(dsep_df_list[[i]])) dsep_df_list[[i]]$Rhat <- round(as.numeric(dsep_df_list[[i]]$Rhat), 3)
        if ("n.eff" %in% names(dsep_df_list[[i]])) dsep_df_list[[i]]$n.eff <- round(as.numeric(dsep_df_list[[i]]$n.eff), 0)
      }
    }
  }
  dsep_df <- if (length(dsep_df_list) > 0) do.call(rbind, dsep_df_list) else NULL
  result <- list(
    dsep = if (!is.null(dsep_df)) list(results = dsep_df) else NULL,
    dsep_tests = current_tests,
    equations = equations,
    data = data,
    original_data = original_data,
    family = if (!is.null(family)) as.list(family) else NULL,
    categorical_vars = attr(data, "categorical_vars"),
    poly_terms = all_poly_terms
  )
  result$call <- original_call
  class(result) <- "because"
  return(result)
}
