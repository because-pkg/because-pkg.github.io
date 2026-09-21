#' Run d-separation tests using JAGS or NIMBLE backends
#'
#' Evaluates the basis set of conditional independence tests using JAGS or NIMBLE,
#' supporting sequential or parallel execution, cross-scale tests, model reuse,
#' and combining resulting MCMC samples.
#'
#' @keywords internal
run_dsep_jags_nimble <- function(
  engine, equations, data, family, structures, structure_obj,
  hierarchical_info, is_hierarchical, random_terms, levels, multiscale,
  link_vars, latent, latent_method, id_col, variability, all_poly_terms,
  fixed_equations_temp, induced_cors_in = NULL, dsep_max_obs, aggregate_crossscale,
  parallel, n.cores, n.chains, n.iter, n.burnin, n.thin, n.adapt,
  ic_recompile, fix_residual_variance, quiet, priors, monitor_mode,
  expand_ordered, nimble_samplers, adapt_delta, max_treedepth,
  prior_scale_fixed, verbose,
  family_obj = NULL, hierarchy = NULL, original_call = NULL,
  original_data = NULL, random = NULL, response_vars_with_na = NULL,
  reuse_models = FALSE, structure = NULL,
  WAIC = FALSE, DIC = FALSE
) {
  monitor <- monitor_mode
  induced_cors <- induced_cors_in

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
    all_vars_in_model <- unique(unlist(lapply(equations, all.vars)))
    extra_random <- parse_global_random(
      random,
      equations,
      all_vars = all_vars_in_model
    )

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

  # Deduplicate tests in case translation created duplicates
  dsep_test_strs <- sapply(dsep_tests, function(eq) {
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
    reused_results <- incremental_check$found
    tests_to_run_indices <- incremental_check$missing_indices

    if (length(tests_to_run_indices) < 2) {
      use_parallel <- FALSE
    }
  }

  if (!quiet) {
    if (use_parallel) {
      message(sprintf(
        "Running %d d-separation tests in parallel on %d cores...",
        length(tests_to_run_indices),
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

  new_results_list <- vector("list", length(dsep_tests))

  if (length(reused_results) > 0) {
    for (i in seq_along(reused_results)) {
      if (!is.null(reused_results[[i]])) {
        new_results_list[[i]] <- reused_results[[i]]
      }
    }
  }

  if (length(tests_to_run_indices) > 0) {
    if (use_parallel) {
      cl <- NULL
      if (is.null(cl)) {
        cl <- parallel::makeCluster(n.cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      }

      parent_exts <- grep("^because\\.", loadedNamespaces(), value = TRUE)

      parallel::clusterEvalQ(cl, {
        if (requireNamespace("because", quietly = TRUE)) {
          library(because)
        }
      })

      if (length(parent_exts) > 0) {
        parallel::clusterExport(cl, "parent_exts", envir = environment())
        parallel::clusterEvalQ(cl, {
          for (ext in parent_exts) {
            if (requireNamespace(ext, quietly = TRUE)) {
              library(ext, character.only = TRUE)
            }
          }
        })
      }

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

      if (requireNamespace("pbapply", quietly = TRUE)) {
        par_results <- pbapply::pblapply(
          X = tests_to_run_indices,
          cl = cl,
          FUN = function(i) {
            test_eq <- dsep_tests[[i]]
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
                variability = variability,
                ic_recompile = ic_recompile,
                dsep_max_obs = dsep_max_obs,
                aggregate_crossscale = aggregate_crossscale
              )
            }, error = function(e) {
              list(
                test_index = i,
                error = e$message,
                test_eq = test_eq
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
                variability = variability,
                ic_recompile = ic_recompile,
                dsep_max_obs = dsep_max_obs,
                aggregate_crossscale = aggregate_crossscale
              )
            }, error = function(e) {
              list(
                test_index = i,
                error = e$message,
                test_eq = test_eq
              )
            })
          }
        )
      }

      for (k in seq_along(tests_to_run_indices)) {
        i <- tests_to_run_indices[k]
        new_results_list[[i]] <- par_results[[k]]
      }
    } else {
      # Sequential execution
      iterator <- tests_to_run_indices
      if (requireNamespace("pbapply", quietly = TRUE) && length(iterator) > 0) {
        seq_results <- pbapply::pblapply(iterator, function(i) {
          test_eq <- dsep_tests[[i]]
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
              variability = variability,
              ic_recompile = ic_recompile,
              dsep_max_obs = dsep_max_obs,
              aggregate_crossscale = aggregate_crossscale
            )
          }, error = function(e) {
            list(
              test_index = i,
              error = e$message,
              test_eq = test_eq
            )
          })
        })

        for (k in seq_along(tests_to_run_indices)) {
          i <- tests_to_run_indices[k]
          new_results_list[[i]] <- seq_results[[k]]
        }
      } else {
        for (i in tests_to_run_indices) {
          test_eq <- dsep_tests[[i]]
          current_monitor <- "interpretable"

          res <- tryCatch({
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
              variability = variability,
              ic_recompile = ic_recompile,
              dsep_max_obs = dsep_max_obs,
              aggregate_crossscale = aggregate_crossscale
            )
          }, error = function(e) {
            list(
              test_index = i,
              error = e$message,
              test_eq = test_eq
            )
          })
          new_results_list[[i]] <- res
        }
      }
    }
  }

  results <- new_results_list
  combined_models <- list()

  for (i in seq_along(results)) {
    res_item <- results[[i]]
    if (is.null(res_item)) next

    if (!is.null(res_item$error)) {
      if (!quiet) {
        warning(paste("D-sep test failed for", deparse(res_item$test_eq), ":", res_item$error))
      }
      next
    }

    tryCatch(
      {
        samples <- res_item$samples
        param_map <- res_item$param_map
        model <- res_item$model

        combined_models[[i]] <- model

        for (ch in 1:coda::nchain(samples)) {
          chain <- samples[[ch]]
          colnames(chain) <- paste0(colnames(chain), "_", i)
          samples[[ch]] <- chain
        }

        if (!is.null(param_map) && nrow(param_map) > 0) {
          param_map$parameter <- paste0(param_map$parameter, "_", i)
        }

        if (is.null(combined_samples)) {
          combined_samples <- samples
        } else {
          if (coda::niter(combined_samples) != coda::niter(samples)) {
            stop("MCMC iteration mismatch between d-sep tests")
          }
          new_samples <- coda::mcmc.list()
          for (ch in 1:coda::nchain(combined_samples)) {
            mat1 <- combined_samples[[ch]]
            mat2 <- samples[[ch]]
            new_mat <- cbind(mat1, mat2)
            new_samples[[ch]] <- coda::mcmc(
              new_mat,
              start = stats::start(mat1),
              thin = coda::thin(mat1)
            )
          }
          combined_samples <- new_samples
        }

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
    dsep_results = results,
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
