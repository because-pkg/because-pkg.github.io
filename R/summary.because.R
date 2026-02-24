#' Summary for Because Model
#'
#' Summarizes the output of a Because model run.
#'
#' @param show_internal Logical. If \code{TRUE}, shows internal parameters created for
#'   deterministically defined nodes (e.g., \code{beta_..._det_...}). Defaults to \code{FALSE}.
#' @param show_nodes Logical. If \code{TRUE}, shows latent node values (e.g., \code{Age[1]}).
#'   Defaults to \code{FALSE} to prevent clutter when \code{monitor="all"}.
#' @param show_random Logical. If \code{TRUE}, shows random effect estimates (e.g., \code{u_...}).
#'   Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[coda]{summary.mcmc}}.
#'
#' @return A summary object containing statistics for the monitored parameters.
#'   If \code{dsep = TRUE} was used in \code{because}, the summary focuses on
#'   the conditional independence tests.
#'
#' @export
summary.because <- function(
    object,
    show_internal = FALSE,
    show_nodes = FALSE,
    show_random = FALSE,
    ...
) {
    # If this was a d-sep run, we want to format the output specifically
    if (!is.null(object$dsep) && object$dsep) {
        tests <- object$dsep_tests
        map <- object$parameter_map

        # Create a results table
        results <- data.frame(
            Test = character(),
            Parameter = character(),
            Estimate = numeric(),
            LowerCI = numeric(),
            UpperCI = numeric(),
            Indep = character(),
            P = numeric(),
            Rhat = numeric(),
            n.eff = numeric(),
            stringsAsFactors = FALSE
        )

        dsep_results <- object$dsep_results

        # Check if results are available
        if (is.null(dsep_results)) {
            warning("No d-separation test results found in model object.")
            return(invisible(NULL))
        }

        # Use lapply for linear performance instead of iterative rbind
        results_list <- lapply(seq_along(dsep_results), function(i) {
            res_i <- dsep_results[[i]]

            # Skip if result is available
            if (is.null(res_i)) {
                return(NULL)
            }

            test_formula <- tests[[i]]
            test_var <- attr(test_formula, "test_var")
            # If test_var is missing from attribute, infer from formula
            if (is.null(test_var)) {
                rhs <- labels(stats::terms(test_formula))
                if (length(rhs) > 0) test_var <- rhs[1]
            }

            response <- as.character(test_formula)[2]

            # Use local map from this specific test run
            map_i <- res_i$param_map

            # Find the parameter name for this path (response ~ test_var)
            param_row <- map_i[
                map_i$response == response &
                    map_i$predictor == test_var,
            ]

            if (nrow(param_row) == 0) {
                # Fallback: scan for any parameter associated with test_var
                param_row <- map_i[map_i$predictor == test_var, ]
            }

            if (nrow(param_row) == 0) {
                warning(paste(
                    "Could not find parameter for test:",
                    format(test_formula)
                ))
                return(NULL)
            }

            param_name <- param_row$parameter[1] # Take first match

            # Optimize: Subset samples to ONLY the relevant parameter before summarizing
            # This is CRITICAL for performance. Otherwise we summarize 1000s of latent nodes per test.
            # We use distinct subsetting that works for mcmc.list

            # Check if parameter exists in samples
            if (!param_name %in% coda::varnames(res_i$samples)) {
                warning(paste(
                    "Parameter",
                    param_name,
                    "not found in samples for test",
                    i
                ))
                return(NULL)
            }

            # Subset to just the parameter of interest
            # Note: mcmc.list subsetting with [, vars] works in standard coda
            sub_samples <- res_i$samples[, param_name, drop = FALSE]

            # Summarize only this parameter
            summ_i <- summary(sub_samples)

            # Handle edge case: single parameter returns vector, not matrix
            if (!is.matrix(summ_i$statistics)) {
                # Convert to matrix format
                param_col_name <- colnames(sub_samples[[1]])[1]
                summ_i$statistics <- matrix(
                    summ_i$statistics,
                    nrow = 1,
                    dimnames = list(param_col_name, names(summ_i$statistics))
                )
                summ_i$quantiles <- matrix(
                    summ_i$quantiles,
                    nrow = 1,
                    dimnames = list(param_col_name, names(summ_i$quantiles))
                )
            }

            # Extract statistics
            if (param_name %in% rownames(summ_i$quantiles)) {
                est <- summ_i$statistics[param_name, "Mean"]
                lower <- summ_i$quantiles[param_name, "2.5%"]
                upper <- summ_i$quantiles[param_name, "97.5%"]

                # Diagnostics (locally computed for this chain list)
                # n.chains for this specific run
                n_chains_i <- coda::nchain(sub_samples)

                # Rhat
                p_rhat <- NA
                if (n_chains_i > 1) {
                    # Handle potential single-parameter Rhat issues
                    p_rhat <- tryCatch(
                        {
                            coda::gelman.diag(
                                sub_samples,
                                multivariate = FALSE
                            )$psrf[param_name, 1]
                        },
                        error = function(e) NA
                    )
                }

                # Effective Size
                p_neff <- NA
                eff_size_i <- tryCatch(
                    {
                        coda::effectiveSize(sub_samples)
                    },
                    error = function(e) NULL
                )

                if (!is.null(eff_size_i) && param_name %in% names(eff_size_i)) {
                    p_neff <- eff_size_i[param_name]
                }

                # Independence check
                indep <- if (lower > 0 || upper < 0) "No" else "Yes"

                # P(~0)
                # Use sub_samples which already contains only the param of interest
                samples_matrix <- as.matrix(sub_samples)
                n_samples <- nrow(samples_matrix)
                param_samples <- samples_matrix[, param_name]
                n_above <- sum(param_samples > 0)
                n_below <- sum(param_samples < 0)
                p_approx_0 <- 2 * min(n_above, n_below) / n_samples

                # Construct label - extract all RHS terms including random effects
                # Use deparse to get full formula, then extract RHS
                formula_str <- paste(deparse(test_formula), collapse = " ")
                # Split on ~ and get RHS
                rhs_full <- sub("^[^~]+~\\s*", "", formula_str)
                # Split on + to get individual terms
                all_terms <- trimws(strsplit(rhs_full, "\\+")[[1]])

                # Separate test_var from conditioning vars
                # Remove test_var from all_terms to get conditioning set
                cond_terms <- all_terms[all_terms != test_var]

                test_str <- paste0(
                    response,
                    " _||_ ",
                    test_var,
                    if (length(cond_terms) > 0) {
                        paste0(" | {", paste(cond_terms, collapse = ","), "}")
                    } else {
                        " | {} "
                    }
                )

                return(data.frame(
                    Test = test_str,
                    Parameter = param_name,
                    Estimate = round(est, 3),
                    LowerCI = round(lower, 3),
                    UpperCI = round(upper, 3),
                    Indep = indep,
                    P = round(p_approx_0, 3), # Renamed from P_approx_0
                    Rhat = round(p_rhat, 3),
                    n.eff = round(p_neff, 0),
                    stringsAsFactors = FALSE
                ))
            }
            return(NULL)
        })

        # Bind all results efficiently
        results <- do.call(rbind, results_list)

        # Store components in list instead of printing
        out <- list(
            type = "dsep",
            results = results
        )
        class(out) <- "summary.because"
        return(out)
    } else {
        # Standard summary
        # Use stored summary if available, otherwise calculate it
        if (!is.null(object$summary)) {
            summ <- object$summary
        } else {
            summ <- summary(object$samples, ...)
        }

        # Calculate convergence diagnostics
        n_chains <- coda::nchain(object$samples)

        # Effective sample size
        eff_size <- tryCatch(
            coda::effectiveSize(object$samples),
            error = function(e) return(NULL)
        )

        # Gelman-Rubin diagnostic (Rhat)
        # Use stored Rhat if available (calculated in because)
        rhat <- NULL
        if ("Rhat" %in% colnames(summ$statistics)) {
            rhat <- summ$statistics[, "Rhat"]
            names(rhat) <- rownames(summ$statistics)
        } else if (n_chains > 1) {
            # Calculate if not stored
            rhat <- tryCatch(
                coda::gelman.diag(object$samples, multivariate = FALSE)$psrf[,
                    1
                ],
                error = function(e) return(NULL)
            )
        }

        # Combine statistics with diagnostics
        stats_table <- summ$statistics[,
            c("Mean", "SD", "Naive SE", "Time-series SE"),
            drop = FALSE
        ]
        quant_table <- summ$quantiles[, c("2.5%", "50%", "97.5%"), drop = FALSE]

        # Create combined table
        combined <- cbind(stats_table, quant_table)

        # Smart Filtering: Hide internal deterministic parameters by default
        if (!show_internal) {
            # Regex to identify internal deterministic link parameters
            # Matches "beta_..._det_..." AND "beta_..._1_times_..." (legacy)
            # The 'beta' is strict because we might want to see 'sigma_det' if that ever exists?
            # Usually these are betas linking the deterministic term.

            # Pattern: ^beta_.*_det_ OR ^beta_.*_times_ (for the pre-fix version)
            # Safe bet: contains "_det_" or "_times_" combined with beta
            internal_idx <- grep("^beta_.*(_det_|_times_)", rownames(combined))

            if (length(internal_idx) > 0) {
                combined <- combined[-internal_idx, , drop = FALSE]
            }
        }

        # Filter Nodes (e.g. Age[1], Rank[5]...)
        # Heuristic: Parameters that overlap with variable names in the data/model, usually indexed.
        # We want to keep 'alpha...', 'beta...', 'sigma...', 'tau...', 'lambda...', 'cutpoint...'
        # We want to hide 'Age[...]', 'Rank[...]' if they are monitored nodes.
        if (!show_nodes) {
            # Keep standard structural prefixes
            # We filter OUT things that look like nodes if they are not structural.

            is_structural <- grepl(
                "^(alpha|beta|sigma|tau|lambda|rho|psi|p_|r_|z_|cutpoint)",
                rownames(combined)
            )
            is_random <- grepl("^u_", rownames(combined))

            # If show_nodes=FALSE, we hide things that are NOT structural and NOT random (unless show_random is on)
            to_keep <- is_structural

            if (show_random) {
                to_keep <- to_keep | is_random
            }

            # Note: Internal parameters (beta_..._det_...) were filtered earlier if show_internal=FALSE.
            # However, if show_internal=TRUE, they start with 'beta', so is_structural keeps them. Perfect.

            combined <- combined[to_keep, , drop = FALSE]
        }

        # Add Rhat if available
        if (!is.null(rhat)) {
            # Ensure alignment
            rhat_aligned <- rhat[rownames(combined)]
            combined <- cbind(combined, Rhat = rhat_aligned)
        }

        # Add n.eff if available
        if (!is.null(eff_size)) {
            eff_aligned <- eff_size[rownames(combined)]
            combined <- cbind(combined, n.eff = eff_aligned)
        }

        # Round output for cleaner display
        # Mean, SD, SE, quantiles: 3 decimal places
        # Rhat: 3 decimal places
        # n.eff: 0 decimal places (integer)
        cols_to_round_3 <- intersect(
            colnames(combined),
            c(
                "Mean",
                "SD",
                "Naive SE",
                "Time-series SE",
                "2.5%",
                "50%",
                "97.5%",
                "Rhat"
            )
        )
        for (col in cols_to_round_3) {
            combined[, col] <- round(combined[, col], 3)
        }
        if ("n.eff" %in% colnames(combined)) {
            combined[, "n.eff"] <- round(combined[, "n.eff"], 0)
        }

        # Store components in list instead of printing
        out <- list(
            type = "standard",
            results = combined,
            DIC = object$DIC,
            WAIC = object$WAIC
        )
        class(out) <- "summary.because"
        return(out)
    }
}

#' Print Summary for Because Model
#'
#' @param x A summary object of class \code{"summary.because"}.
#' @param ... Additional arguments.
#'
#' @export
print.summary.because <- function(x, ...) {
    if (x$type == "dsep") {
        results <- x$results
        cat("d-separation Tests\n")
        cat("==================\n\n")
        # Print each test on a separate block
        for (i in 1:nrow(results)) {
            cat(paste0("Test: ", results$Test[i]), "\n")
            # Print stats row without the Test column
            print(results[i, -1], row.names = FALSE)
            cat("\n")
        }

        cat("\nLegend:\n")
        cat(
            "  Indep: 'Yes' = Conditionally Independent, 'No' = Dependent (based on 95% CI)\n"
        )
        cat(
            "  P: Bayesian probability that the posterior distribution overlaps with zero\n"
        )
    } else {
        # Standard summary
        print(x$results)

        if (!is.null(x$DIC)) {
            cat("\nDIC:\n")
            print(x$DIC)
        }

        if (!is.null(x$WAIC)) {
            cat("\nWAIC:\n")
            print(x$WAIC)
        }
    }
    invisible(x)
}
