#' Compute MCMC Summary and Gelman-Rubin R-hat Statistic
#'
#' @param samples A coda mcmc.list object.
#' @param n.chains Number of MCMC chains.
#' @return A summary object with R-hat appended to statistics if n.chains > 1,
#'   excluding internal log_lik parameters.
#' @noRd
compute_mcmc_summary <- function(samples, n.chains) {
  if (is.null(samples)) {
    return(NULL)
  }

  sum_stats <- summary(samples)

  # Explicitly calculate R-hat if multiple chains
  if (n.chains > 1) {
    tryCatch(
      {
        # Manual R-hat calculation to avoid coda::gelman.diag issues
        # with parallel chains and R scoping problems
        n_chains <- length(samples)

        # Use base R colnames to get parameter names
        first_chain <- as.matrix(samples[[1]])
        pnames <- colnames(first_chain)
        n_params <- length(pnames)

        psrf <- matrix(NA, nrow = n_params, ncol = 2)
        rownames(psrf) <- pnames
        colnames(psrf) <- c("Point est.", "Upper C.I.")

        # Convert chains to matrices ONCE outside the loop
        chain_matrices <- lapply(samples, as.matrix)

        for (idx in seq_len(n_params)) {
          # Extract column idx from each chain matrix
          vals <- do.call(cbind, lapply(chain_matrices, function(m) m[, idx]))

          # Check for constant chains (variance 0)
          # Use explicit variance calculation to avoid R scoping issues
          chain_vars <- numeric(ncol(vals))
          for (col_idx in seq_len(ncol(vals))) {
            chain_vars[col_idx] <- stats::var(vals[, col_idx])
          }

          if (any(chain_vars < 1e-10)) {
            psrf[idx, 1] <- 1.0 # If constant, Rhat is 1
            next
          }

          # Calculate B/W (Gelman-Rubin statistic)
          n_samples <- nrow(vals)
          chain_means <- colMeans(vals)

          # Between-chain variance
          B <- n_samples * stats::var(chain_means)

          # Within-chain variance
          W <- mean(chain_vars)

          # Estimated variance
          var_plus <- (n_samples - 1) / n_samples * W + B / n_samples

          # R-hat
          rhat <- sqrt(var_plus / W)
          psrf[idx, 1] <- rhat
        }

        # Match parameter names
        common_params <- intersect(
          rownames(sum_stats$statistics),
          rownames(psrf)
        )

        if (length(common_params) > 0) {
          sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = NA)
          rhat_col_idx <- which(colnames(sum_stats$statistics) == "Rhat")

          # Use numeric indexing to avoid any strange symbol evaluation
          for (j in seq_along(common_params)) {
            p <- common_params[j]
            row_idx <- which(rownames(sum_stats$statistics) == p)
            psrf_row_idx <- which(rownames(psrf) == p)
            if (length(row_idx) == 1 && length(psrf_row_idx) == 1) {
              sum_stats$statistics[row_idx, rhat_col_idx] <- psrf[
                psrf_row_idx,
                1
              ]
            }
          }
        }
      },
      error = function(e) {
        warning("Could not calculate R-hat: ", e$message)
      }
    )
  }

  # Filter internal parameters (log_lik) from summary parameters
  # We keep them in samples for WAIC calculation but hide them from the summary output
  if (!is.null(sum_stats)) {
    if (is.matrix(sum_stats$statistics)) {
      rows_to_keep <- !grepl("^log_lik", rownames(sum_stats$statistics))
      sum_stats$statistics <- sum_stats$statistics[
        rows_to_keep,
        ,
        drop = FALSE
      ]
      sum_stats$quantiles <- sum_stats$quantiles[rows_to_keep, , drop = FALSE]
    } else {
      # Single parameter case (statistics is a vector)
      # Check if the single parameter is log_lik
      param_name <- colnames(samples[[1]])
      if (length(param_name) == 1 && grepl("^log_lik", param_name)) {
        sum_stats <- NULL
      }
    }
  }

  sum_stats
}


#' Assemble Final because S3 Object
#'
#' Packages model, samples, summaries, information criteria (DIC/WAIC), and metadata.
#'
#' @noRd
assemble_because_result <- function(
  model,
  model_code,
  model_file,
  samples,
  data,
  original_data,
  equations,
  random,
  random_terms,
  structure,
  structures,
  latent,
  distribution,
  family,
  variability,
  all_poly_terms,
  dsep,
  dsep_tests,
  dsep_results,
  parameter_map,
  induced_cors,
  scale_info,
  stack_res = NULL,
  engine = "jags",
  saved_nimble_compiled = NULL,
  saved_nimble_cmodel = NULL,
  saved_nimble_samplers = NULL,
  nimble_waic = NULL,
  parallel = FALSE,
  n.cores = 1,
  n.chains = 3,
  n.iter = 2000,
  n.burnin = 1000,
  n.thin = 1,
  n.adapt = 1000,
  DIC = FALSE,
  WAIC = FALSE,
  ic_recompile = FALSE,
  extension_inits = list(),
  quiet = FALSE,
  id_col = NULL,
  original_call = NULL,
  hierarchical_info = NULL,
  monitor = NULL
) {
  # Summarize posterior
  sum_stats <- compute_mcmc_summary(samples, n.chains)

  # Initialize result object
  result <- list(
    model = model,
    model_code = model_code,
    data = data, # Store data for recompilation if needed
    input = list(
      equations = equations,
      random = random,
      structure = structure,
      data = original_data, # Store original data too for safety
      latent = latent,
      distribution = distribution,
      family = if (!is.null(family)) as.list(family) else NULL,
      variability = variability,
      poly_terms = all_poly_terms # Needed by plot_dag to reconstruct diamond nodes
    ),
    samples = samples,
    summary = sum_stats,
    monitor = monitor,
    modfile = model_file,
    dsep = dsep,
    dsep_tests = dsep_tests,
    dsep_results = dsep_results,
    parameter_map = parameter_map,
    induced_correlations = induced_cors,
    scale_info = scale_info,
    stacked_data = if (!is.null(stack_res) && isTRUE(stack_res$is_stacked)) {
      stack_res$data
    } else {
      NULL
    }
  )

  result$engine <- engine
  if (engine == "nimble") {
    result$nimble_compiled <- saved_nimble_compiled
    result$nimble_cmodel   <- saved_nimble_cmodel
    result$nimble_samplers <- saved_nimble_samplers
  }
  result$samples        <- samples
  result$parameter_map  <- parameter_map
  result$data           <- data
  result$original_data  <- original_data
  result$family         <- family
  result$categorical_vars <- attr(data, "categorical_vars")
  result$poly_terms     <- all_poly_terms
  result$equations      <- equations
  result$parallel       <- parallel
  result$n.cores        <- n.cores

  # --- Result Enrichment ---
  # If we have a structure, try to extract labels for ordering
  if (!is.null(structure)) {
    result$species_order <- get_order_labels_hook(structure)
  } else if (!is.null(id_col) && is.data.frame(original_data)) {
    # If no structure but ID col provided
    result$species_order <- as.character(original_data[[id_col]])
  }

  # Assign class immediately (needed for print/summary/waic methods)
  result$call <- original_call
  class(result) <- "because"

  # Add DIC and WAIC
  # For parallel runs, recompile the model if ic_recompile=TRUE
  if (
    (DIC || WAIC) && parallel && n.cores > 1 && n.chains > 1 && ic_recompile
  ) {
    message("Recompiling model for DIC/WAIC calculation...")

    # Recompile model with 2 chains for IC calculation (DIC requires >=2)
    ic_inits <- lapply(1:2, function(i) {
      c(
        extension_inits,
        list(
          .RNG.name = "base::Wichmann-Hill",
          .RNG.seed = 12345 + i
        )
      )
    })

    ic_model <- rjags::jags.model(
      model_file,
      data = data,
      inits = ic_inits,
      n.chains = 2,
      n.adapt = n.adapt,
      quiet = quiet
    )

    # Short burn-in (use a fraction of original)
    if (n.burnin > 0) {
      update(ic_model, n.iter = min(n.burnin, 500))
    }

    # Compute DIC
    if (DIC) {
      if (n.iter > n.burnin) {
        result$DIC <- rjags::dic.samples(
          ic_model,
          n.iter = min(n.iter - n.burnin, 1000)
        )
      } else {
        result$DIC <- NULL
      }
    }
  } else if ((DIC || WAIC) && parallel && n.cores > 1 && n.chains > 1) {
    # Parallel without recompilation - warn user
    if (DIC) {
      warning(
        "DIC calculation disabled for parallel chains. Set ic_recompile=TRUE to compute DIC."
      )
      result$DIC <- NULL
    }
  } else {
    # Sequential execution - use standard approach
    if (DIC) {
      if (engine == "jags") {
        if (n.iter > n.burnin) {
          result$DIC <- rjags::dic.samples(model, n.iter = n.iter - n.burnin)
        } else {
          result$DIC <- NULL
        }
      } else {
        result$DIC <- NULL # NIMBLE does not use rjags::dic.samples
      }
    }
  }

  # Compute WAIC if requested (must be after class assignment)
  if (WAIC) {
    if (engine == "nimble" && !is.null(nimble_waic)) {
      # NIMBLE built-in WAIC: lppd = log pointwise predictive density (no penalty),
      # pWAIC = effective parameters, WAIC = -2*(lppd - pWAIC).
      # elpd_waic (as in because) = lppd - pWAIC  (NOT just lppd)
      nimble_elpd     <- nimble_waic$lppd - nimble_waic$pWAIC
      nimble_pwaic    <- nimble_waic$pWAIC
      nimble_waic_val <- nimble_waic$WAIC   # = -2 * nimble_elpd
      # n_obs: N observations x number of modelled response variables
      n_obs_nimble <- tryCatch({
        n_resp <- length(unique(result$parameter_map$response))
        as.integer(data$N * n_resp)
      }, error = function(e) NA_integer_)
      n_samples_nimble <- as.integer((n.iter - n.burnin) / n.thin) * n.chains
      waic_df <- data.frame(
        Estimate = c(nimble_elpd, nimble_pwaic, nimble_waic_val),
        SE = c(NA_real_, NA_real_, NA_real_), # no pointwise SE from NIMBLE built-in WAIC
        row.names = c("elpd_waic", "p_waic", "waic")
      )
      attr(waic_df, "dims") <- c(n_obs = n_obs_nimble, n_samples = n_samples_nimble)
      class(waic_df) <- c("because_waic", "data.frame")
      result$WAIC <- waic_df
    } else {
      result$WAIC <- because_waic(result)
    }
  }

  # Preserve hierarchical metadata for diagnostics even if data was flat
  result$hierarchical_info <- hierarchical_info

  return(result)
}
