#' Run d-separation tests for a because model
#'
#' This function contains the full d-separation testing pipeline extracted
#' from because(). It handles both JAGS/NIMBLE and NumPyro engines, parallel
#' execution, cross-scale dispatch, and latent variable methods.
#'
#' All arguments match the corresponding local variables inside because().
#' @keywords internal
run_because_dsep <- function(
  dsep, engine, equations, data, family, structures, structure_obj,
  hierarchical_info, is_hierarchical, random_terms, levels, multiscale,
  link_vars, latent, latent_method, id_col, variability, all_poly_terms,
  fixed_equations_temp, induced_cors_in, dsep_max_obs, aggregate_crossscale,
  parallel, n.cores, n.chains, n.iter, n.burnin, n.thin, n.adapt,
  ic_recompile, fix_residual_variance, quiet, priors, monitor_mode,
  expand_ordered, nimble_samplers, adapt_delta, max_treedepth,
  prior_scale_fixed, verbose
) {
  dsep_tests <- NULL
  induced_cors <- induced_cors_in


  if (dsep) {
    if (engine == "numpyro") {
      eq_strings <- sapply(equations, function(eq) paste(deparse(eq), collapse=" "))
      
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
                mean = round(mean(vec, na.rm=TRUE), 3),
                ci_2.5 = round(unname(quantile(vec, 0.025, na.rm=TRUE)), 3),
                ci_97.5 = round(unname(quantile(vec, 0.975, na.rm=TRUE)), 3),
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
    # Extension Hook: Remove specialized variables from potential latents

    # Force WAIC and DIC off for d-separation testing (not needed for conditional independence tests)
    if (WAIC || DIC) {
      if (!quiet) {
        message(
          "Note: WAIC and DIC are not computed for d-separation tests (not needed for conditional independence testing)."
        )
      }
      WAIC <- FALSE
      DIC <- FALSE
    }


    # Expanded random terms for d-sep (include all variables as potential responses)
    # This ensures that root nodes in the DAG get random effects in d-sep tests
    random_terms_for_dsep <- random_terms
    if (!is.null(random)) {
      # Use all variables currently in the model equations
      all_vars_in_model <- unique(unlist(lapply(equations, all.vars)))
      extra_random <- parse_global_random(
        random,
        equations,
        all_vars = all_vars_in_model
      )

      # Combine and deduplicate
      combined_rand <- c(random_terms, extra_random)
      rand_keys <- sapply(combined_rand, function(x) {
        paste(x$response, x$group, sep = "|")
      })
      random_terms_for_dsep <- combined_rand[!duplicated(rand_keys)]
    }

    dsep_result <- because_dsep(
      equations,
      latent = latent,
      random_terms = random_terms_for_dsep,
      categorical_vars = if (!is.null(attr(data, "categorical_vars"))) {
        attr(data, "categorical_vars")
      } else {
        NULL
      },
      family = if (!is.null(family)) as.list(family) else NULL,
      poly_terms = all_poly_terms,
      quiet = quiet,
      hierarchical_info = hierarchical_info
    )

    # Extract tests and correlations
    if (!is.null(latent)) {
      dsep_tests <- dsep_result$tests
      induced_cors <- dsep_result$correlations

      if (!quiet && length(induced_cors) > 0) {
        message(
          "Found ",
          length(induced_cors),
          " induced correlation(s) from latent variable(s)"
        )
      }
    } else {
      dsep_tests <- dsep_result
    }

    # Extension Hook: Translate tests where response is specialized (e.g. psi_)
    dsep_tests <- lapply(dsep_tests, function(eq) {
      attr_val <- attr(eq, "test_var")
      new_eq <- dsep_test_translation_hook(family_obj, eq)
      if (!is.null(attr_val)) attr(new_eq, "test_var") <- attr_val
      new_eq
    })

    # (Original translation block removed)

    # Deduplicate tests in case translation created duplicates
    dsep_test_strs <- sapply(dsep_tests, function(eq) {
      # Use trimws and a single space to normalize for string comparison
      str <- paste(deparse(eq), collapse = " ")
      gsub("\\s+", " ", trimws(str))
    })

    dsep_tests <- dsep_tests[!duplicated(dsep_test_strs)]

    if (length(dsep_tests) == 0) {
      stop(
        "No d-separation tests implied by the model (model is saturated). Stopping run."
      )
    }

    # --- A priori cross-hierarchy filter ---
    # Tests where the response and the focal predictor belong to orthogonal
    # hierarchical branches (e.g., species-level vs. site-level) are trivially
    # satisfied by construction: there is no causal path between the two branches
    # except through the observation level, and no cross-level index exists in
    # the data to estimate such a regression. Skip them with a clear message
    # rather than letting JAGS fail with "Unknown variable species_idx_site".
    cross_hier_flags <- sapply(dsep_tests, function(eq) {
      is_cross_hierarchy_test(eq, hierarchical_info)
    })
    n_cross <- sum(cross_hier_flags)
    if (n_cross > 0 && !quiet) {
      skipped_labels <- sapply(dsep_tests[cross_hier_flags], function(eq) {
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
    dsep_tests_skipped <- dsep_tests[cross_hier_flags]
    dsep_tests <- dsep_tests[!cross_hier_flags]

    if (length(dsep_tests) == 0) {
      warning("All d-separation tests were cross-hierarchy (trivially satisfied). ",
              "No JAGS models will be run.")
    }

    # Run tests sequentially to avoid cyclic dependencies in JAGS
    # Decide whether to run tests in parallel
    use_parallel <- parallel && n.cores > 1 && length(dsep_tests) > 1

    # INCREMENTAL D-SEP: Check for reusable tests
    reused_results <- list()
    tests_to_run_indices <- seq_along(dsep_tests)

    if (!is.null(reuse_models)) {
      incremental_check <- find_reusable_tests(
        dsep_tests,
        equations,
        reuse_models,
        data,
        family = if (!is.null(family)) as.list(family) else NULL,
        quiet = quiet
      )
      reused_results <- incremental_check$found # List with results at matching indices, NULL otherwise
      tests_to_run_indices <- incremental_check$missing_indices

      # Override parallel decision if few tests remain
      if (length(tests_to_run_indices) < 2) {
        use_parallel <- FALSE
      }
    }

    if (!quiet) {
      if (use_parallel) {
        message(sprintf(
          "Running %d d-separation tests in parallel on %d cores...",
          length(tests_to_run_indices), # Corrected message to show actual runs
          n.cores
        ))
      } else {
        message(sprintf(
          "Running %d d-separation tests sequentially...",
          length(tests_to_run_indices)
        ))
      }
    }

    combined_samples <- NULL
    combined_map <- NULL

    # Define function to run a single d-sep test
    # run_single_dsep_test_v2 removed from here (now top-level)

    # Run tests (parallel or sequential)
    # Only loop over tests_to_run_indices
    # We still need to populate a full results list of length(dsep_tests)
    # The 'reused_results' list has NULLs for missing tests and values for found ones.

    new_results_list <- vector("list", length(dsep_tests))

    # Fill in reused results first
    if (length(reused_results) > 0) {
      for (i in seq_along(reused_results)) {
        if (!is.null(reused_results[[i]])) {
          new_results_list[[i]] <- reused_results[[i]]
        }
      }
    }

    if (length(tests_to_run_indices) > 0) {
      if (use_parallel) {
        # Setup cluster if not provided
        if (is.null(cl)) {
          cl <- parallel::makeCluster(n.cores)
          on.exit(parallel::stopCluster(cl), add = TRUE)
        }

        # Identify because extensions loaded in the master session
        parent_exts <- grep("^because\\.", loadedNamespaces(), value = TRUE)

        # Ensure workers have the library loaded
        parallel::clusterEvalQ(cl, {
          if (requireNamespace("because", quietly = TRUE)) {
            library(because)
          }
        })

        # Load detected extensions on workers
        if (length(parent_exts) > 0) {
            # Pass the list of extensions to workers
            parallel::clusterExport(cl, "parent_exts", envir = environment())
            parallel::clusterEvalQ(cl, {
                for (ext in parent_exts) {
                    if (requireNamespace(ext, quietly = TRUE)) {
                        library(ext, character.only = TRUE)
                    }
                }
            })
        }

        # Export necessary objects to cluster
        parallel::clusterExport(
          cl,
          c(
            "original_data",
            "structure",
            "monitor",
            "n.chains",
            "n.iter",
            "n.burnin",
            "n.thin",
            "n.adapt",
            "variability",
            "family",
            "latent",
            "latent_method",
            "ic_recompile",
            "equations",
            "random_terms",
            "hierarchical_info",
            "levels",
            "hierarchy",
            "multiscale",
            "link_vars",
            "fix_residual_variance",
            "run_single_dsep_test_v2",
            "dsep_tests",
            "extract_random_effects",
            "engine",
            "nimble_samplers",
            "because",
            "quiet",
            "random",
            "get_data_for_variables",
            "infer_variable_level",
            "get_level_depth",
            "sanitize_term_name",
            "dsep_tree_hook",
            "dsep_equations_hook"
          ),
          envir = environment()
        )

        # Run tests in parallel
        # Only run the necessary tests
        # We use pbapply for a live progress bar if available, otherwise fallback to parLapply
        # Note: pblapply and parLapply have different argument orders for 'cl'
        if (requireNamespace("pbapply", quietly = TRUE)) {
          par_results <- pbapply::pblapply(
            X = tests_to_run_indices,
            cl = cl,
            FUN = function(i) {
            test_eq <- dsep_tests[[i]]

            # Use "interpretable" monitoring mode.
            current_monitor <- "interpretable"

            tryCatch({
                run_single_dsep_test_v2(
                  i,
                  test_eq,
                  current_monitor,
                  engine = engine,
                  nimble_samplers = nimble_samplers,
                  quiet = quiet,
                  original_data = original_data,
                  hierarchical_info = hierarchical_info,
                  random_terms = random_terms,
                  equations = equations,
                family = if (!is.null(family)) as.list(family) else NULL,
                structure = structure,
                levels = levels,
                hierarchy = hierarchy,
                multiscale = multiscale,
                link_vars = link_vars,
                fix_residual_variance = fix_residual_variance,
                latent = latent,
                latent_method = latent_method,
                n.chains = n.chains,
                n.iter = n.iter,
                n.burnin = n.burnin,
                n.thin = n.thin,
                n.adapt = n.adapt,
                ic_recompile = ic_recompile,
                random = random,
                id_col = id_col,
                variability = variability,
                dsep_max_obs = dsep_max_obs,
                aggregate_crossscale = aggregate_crossscale
              )
            }, error = function(e) {
              # Return error info to parent session for reporting
              list(
                error = conditionMessage(e),
                test_idx = i,
                equation = paste(deparse(test_eq), collapse = " ")
              )
            })
          }
        )
        } else {
          par_results <- parallel::parLapply(
            cl = cl,
            X = tests_to_run_indices,
            fun = function(i) {
              test_eq <- dsep_tests[[i]]

              # Use "interpretable" monitoring mode.
              current_monitor <- "interpretable"

              tryCatch({
                run_single_dsep_test_v2(
                  i,
                  test_eq,
                  current_monitor,
                  engine = engine,
                  nimble_samplers = nimble_samplers,
                  quiet = quiet,
                  original_data = original_data,
                  hierarchical_info = hierarchical_info,
                  random_terms = random_terms,
                  equations = equations,
                  family = if (!is.null(family)) as.list(family) else NULL,
                  structure = structure,
                  levels = levels,
                  hierarchy = hierarchy,
                  multiscale = multiscale,
                  link_vars = link_vars,
                  fix_residual_variance = fix_residual_variance,
                  latent = latent,
                  latent_method = latent_method,
                  n.chains = n.chains,
                  n.iter = n.iter,
                  n.burnin = n.burnin,
                  n.thin = n.thin,
                  n.adapt = n.adapt,
                  ic_recompile = ic_recompile,
                  random = random,
                  id_col = id_col,
                  variability = variability,
                  dsep_max_obs = dsep_max_obs,
                  aggregate_crossscale = aggregate_crossscale
                )
              }, error = function(e) {
                # Return error info to parent session for reporting
                list(
                  error = conditionMessage(e),
                  test_idx = i,
                  equation = paste(deparse(test_eq), collapse = " ")
                )
              })
            }
          )
        }

        # Merge parallel results into new_results_list
        for (j in seq_along(tests_to_run_indices)) {
          idx <- tests_to_run_indices[j]
          res <- par_results[[j]]
          
          # Check if the result is an error object from tryCatch
          if (is.list(res) && !is.null(res$error)) {
            if (!quiet) {
              warning(sprintf(
                "D-sep test %d/%d skipped due to error: %s\n  Equation: %s",
                res$test_idx, length(dsep_tests), res$error, res$equation
              ))
            }
            new_results_list[[idx]] <- NULL
          } else {
            new_results_list[[idx]] <- res
          }
        }
      } else {
        # Sequential execution
        # Loop over only the needed indices
        seq_results <- lapply(tests_to_run_indices, function(i) {
          test_eq <- dsep_tests[[i]]
          # Monitor betas for ALL predictors in the d-sep equation
          current_monitor <- "interpretable"

          if (!quiet) {
            message(sprintf(
              "  Test %d/%d: %s",
              i,
              length(dsep_tests),
              deparse(test_eq)
            ))
          }
          new_results_list[[i]] <- tryCatch(
            run_single_dsep_test_v2(
              i,
              test_eq,
              current_monitor,
              engine = engine,
              nimble_samplers = nimble_samplers,
              quiet = quiet,
              original_data = original_data,
              hierarchical_info = hierarchical_info,
              random_terms = random_terms,
              equations = equations,
              family = if (!is.null(family)) as.list(family) else NULL,
              structure = structure,
              levels = levels,
              hierarchy = hierarchy,
              multiscale = multiscale,
              link_vars = link_vars,
              fix_residual_variance = fix_residual_variance,
              latent = latent,
              latent_method = latent_method,
              n.chains = n.chains,
              n.iter = n.iter,
              n.burnin = n.burnin,
              n.thin = n.thin,
              n.adapt = n.adapt,
              ic_recompile = ic_recompile,
              random = random,
              id_col = id_col,
              variability = variability,
              dsep_max_obs = dsep_max_obs,
              aggregate_crossscale = aggregate_crossscale
            ),
            error = function(e) {
              if (!quiet) {
                warning(sprintf(
                  "D-sep test %d/%d skipped due to error: %s\n  Equation: %s",
                  i, length(dsep_tests), conditionMessage(e), deparse(test_eq)
                ))
              }
              NULL  # Return NULL for this test
            }
          )
        })

        # Assign sequential results to the main list
        for (j in seq_along(tests_to_run_indices)) {
          idx <- tests_to_run_indices[j]
          new_results_list[[idx]] <- seq_results[[j]]
        }
      }
    } # End if tests_to_run > 0

    results <- new_results_list

    # Combine results
    combined_models <- list()

    for (res_item in results) {
      if (is.null(res_item)) {
        next
      }

      tryCatch(
        {
          samples <- res_item$samples
          param_map <- res_item$param_map
          model_string <- res_item$model
          i <- res_item$test_index

          # Store model for this test
          combined_models[[i]] <- model_string

          # Rename parameters to include equation index to avoid collisions
          # e.g., betaRS becomes betaRS_1 for equation 1, betaRS_2 for equation 2
          for (ch in seq_along(samples)) {
            chain <- samples[[ch]]
            colnames_orig <- colnames(chain)

            # Add suffix _i to all beta, alpha, lambda, tau, rho parameters
            new_colnames <- sapply(colnames_orig, function(name) {
              if (grepl("^(beta|alpha|lambda|tau|rho|sigma)", name)) {
                paste0(name, "_", i)
              } else {
                name
              }
            })

            colnames(chain) <- new_colnames
            samples[[ch]] <- chain
          }

          # Update parameter_map to reflect new names
          if (!is.null(param_map) && nrow(param_map) > 0) {
            param_map$parameter <- paste0(param_map$parameter, "_", i)
          }

          # Combine samples (cbind chains)
          if (is.null(combined_samples)) {
            combined_samples <- samples
          } else {
            # Check if dimensions match
            if (coda::niter(combined_samples) != coda::niter(samples)) {
              stop("MCMC iteration mismatch between d-sep tests")
            }
            # Combine chains: for each chain, cbind the variables
            new_samples <- coda::mcmc.list()
            for (ch in 1:coda::nchain(combined_samples)) {
              # Combine matrices
              mat1 <- combined_samples[[ch]]
              mat2 <- samples[[ch]]
              # All columns from mat2 should be new (due to renaming)
              new_mat <- cbind(mat1, mat2)
              new_samples[[ch]] <- coda::mcmc(
                new_mat,
                start = stats::start(mat1),
                thin = coda::thin(mat1)
              )
            }
            combined_samples <- new_samples
          }

          # Combine parameter maps
          if (is.null(combined_map)) {
            combined_map <- param_map
          } else {
            combined_map <- rbind(combined_map, param_map)
          }
        },
        error = function(e) {
          if (!quiet) {
            warning(paste(
              "D-sep result combination error for test index",
              res_item$test_index,
              ":",
              e$message
            ))
          }
        }
      )
    }

    # Return combined result
    result <- list(
      samples = combined_samples,
      parameter_map = combined_map,
      models = combined_models,
      dsep = TRUE,
      dsep_tests = dsep_tests,
      dsep_results = results, # Store individual test results for summary
      induced_correlations = induced_cors,
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

  # Handle latent variable method
  if (!is.null(latent)) {
    latent_method <- match.arg(latent_method, c("correlations", "explicit"))

    # Force MAG approach when doing d-separation testing
    if (dsep && latent_method == "explicit") {
      if (!quiet) {
        message(
          "Note: d-separation testing with latent variables requires MAG approach. ",
          "Using latent_method = 'correlations'."
        )
      }
      latent_method <- "correlations"
    }

    if (latent_method == "correlations") {
      # MAG approach: marginalize latents, use induced correlations
      # If not already computed by dsep, compute now
      if (is.null(induced_cors)) {
        dsep_result <- because_dsep(
          equations,
          latent = latent,
          random_terms = random_terms,
          hierarchical_info = hierarchical_info,
          family = if (!is.null(family)) as.list(family) else NULL,
          quiet = !dsep
        )
        induced_cors <- dsep_result$correlations
      }

      # Filter out equations involving latent variables
      # Handle equations involving latent variables
      original_eq_count <- length(equations)
      equations <- lapply(equations, function(eq) {
        vars <- all.vars(eq)
        if (any(vars %in% latent)) {
          # Remove latent variables from the formula
          # We construct a new formula string excluding latent vars
          rhs <- labels(terms(eq))
          keep_terms <- rhs[!rhs %in% latent]

          if (length(keep_terms) == 0) {
            # Becomes intercept-only model
            new_eq <- as.formula(paste(as.character(eq)[2], "~ 1"))
          } else {
            # Keep observed predictors
            new_eq <- as.formula(paste(
              as.character(eq)[2],
              "~",
              paste(keep_terms, collapse = " + ")
            ))
          }
          return(new_eq)
        }
        return(eq)
      })

      # We no longer filter them out, so we don't report "removed equations"
      if (!quiet) {
        message(
          "Using MAG approach: marginalized latent variables from structural equations."
        )
      }

      # Check if variables with induced correlations need intercept-only models
      # This happens when a variable is involved in an induced correlation
      # but is not a response variable in any remaining equation
      if (length(induced_cors) > 0) {
        # Get all variables from induced correlations
        vars_with_correlations <- unique(unlist(induced_cors))

        # Get response variables from remaining equations
        response_vars <- sapply(equations, function(eq) {
          as.character(eq)[2]
        })

        # Find variables that need intercept models
        vars_needing_intercept <- setdiff(
          vars_with_correlations,
          response_vars
        )

        if (length(vars_needing_intercept) > 0) {
          # Create intercept-only models: X ~ 1
          intercept_equations <- lapply(vars_needing_intercept, function(v) {
            as.formula(paste(v, "~ 1"))
          })

          # Add to equations list
          equations <- c(equations, intercept_equations)

          if (!quiet) {
            message(
              "Created intercept-only models for ",
              length(vars_needing_intercept),
              " variable(s) with induced correlations: ",
              paste(vars_needing_intercept, collapse = ", ")
            )
          }
        }
      }

      if (!quiet && length(induced_cors) > 0) {
        message(
          "Estimating ",
          length(induced_cors),
          " induced correlation(s) from latent variable(s)"
        )
      }
    } else {
      # Explicit approach: keep all equations, don't use induced correlations
      induced_cors <- NULL

      if (!quiet) {
        message("Using explicit latent variable modeling")
      }
    }
  }

  # Auto-expand categorical variables in equations
  if (!is.null(attr(data, "categorical_vars"))) {
    categorical_vars <- attr(data, "categorical_vars")
    new_equations <- list()

    # Use for loop instead of lapply to safely modify data and collect equations
    for (idx in seq_along(equations)) {
      eq <- equations[[idx]]
      # Parse formula to get all variables
      vars <- all.vars(eq)

      # Check if any predictors are categorical
      for (var in vars) {
        if (var %in% names(categorical_vars)) {
          # Check if var is the response (LHS)
          lhs_var <- all.vars(eq[[2]])
          if (var %in% lhs_var) {
            # Skip expansion if it's the response
            next
          }

          # Get dummy variable names
          levels <- categorical_vars[[var]]$levels
          dummies <- categorical_vars[[var]]$dummies

          # If the parent variable is being imputed (is in response_vars_with_na),
          # we MUST define the dummies deterministically in JAGS to link them.
          if (var %in% response_vars_with_na) {
            for (k in 2:length(levels)) {
              dummy_name <- paste0(var, "_", levels[k])

              # Create deterministic equation: dummy ~ I(var == k)
              det_eq_str <- sprintf("%s ~ I(%s == %d)", dummy_name, var, k)
              det_eq <- stats::as.formula(det_eq_str)

              # Only add if not already present
              eq_exists <- any(sapply(c(equations, new_equations), function(e) {
                deparse(e) == deparse(det_eq)
              }))

              if (!eq_exists) {
                new_equations <- c(new_equations, list(det_eq))
                # Also remove the dummy from the data list so JAGS uses the definition
                if (dummy_name %in% names(data)) {
                  data[[dummy_name]] <- NULL
                }
              }
            }
          }

          # [FIX] Skip substitution if we are inside a deterministic identity definition
          # for an imputed factor (e.g., Rate_2 ~ I(Rate == 2)).
          # This prevents circular expansion of 'Rate' into its own dummies.
          is_identity_mapping <- length(vars) == 2 && grepl(sprintf("I\\(%s == [0-9]+\\)", var), deparse(eq))
          if (is_identity_mapping) {
            next
          }

          # Convert formula to character for manipulation
          eq_str <- paste(deparse(eq), collapse = " ")

          # Replace categorical variable with its dummies (wrapped in parentheses for proper interaction expansion like A*B)
          pattern <- paste0("\\b", var, "\\b")
          replacement <- paste0("(", paste(dummies, collapse = " + "), ")")
          eq_str <- gsub(pattern, replacement, eq_str)

          # Convert back to formula
          eq <- as.formula(eq_str)
          equations[[idx]] <- eq

          if (!quiet) {
            message(sprintf(
              "Expanded '%s' to: %s",
              var,
              paste(dummies, collapse = ", ")
            ))
          }
        }
      }
    }
    # Add the new deterministic equations
    equations <- c(equations, new_equations)
  }


  return(list(dsep_tests = dsep_tests, induced_cors = induced_cors, equations = equations))
}
