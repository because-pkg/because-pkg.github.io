#' @title Identify Reusable D-Separation Tests
#' @description
#' Scans a list of previously fitted 'because' models to find d-separation tests that have already been run
#' and can be reused. Reusability requires strict data identity and identical test formula + required support equations.
#'
#' @param current_tests List of formula objects representing the d-separation tests to run.
#' @param current_equations List of model formulas for the current model (support equations).
#' @param reuse_models List of 'because' model objects to scan for reusable results.
#' @param current_data The data object (list or environment) being used for the current run.
#' @param family Named vector of families (needed for dependency resolution).
#' @param quiet Logical; suppress messages (default = FALSE).
#'
#' @return A list with two components:
#' \itemize{
#'   \item \code{found}: A list of d-separation result objects found in previous models.
#'   \item \code{missing_indices}: A numeric vector of indices in \code{current_tests} that need to be run.
#' }
#'
#' @keywords internal
find_reusable_tests <- function(
    current_tests,
    current_equations,
    reuse_models,
    current_data,
    family = NULL,
    quiet = FALSE
) {
    if (is.null(reuse_models) || length(reuse_models) == 0) {
        return(list(found = list(), missing_indices = seq_along(current_tests)))
    }

    found <- vector("list", length(current_tests))
    found_indices <- integer(0)
    n_reused <- 0

    if (!quiet) {
        message("Scanning previous models for reusable d-separation tests...")
    }

    # Filter valid models first.
    # JAGS/NIMBLE store individual test results in $dsep_results (list of coda objects).
    # NumPyro stores a pre-formatted data.frame in $dsep$results.
    valid_models <- list()
    for (m_idx in seq_along(reuse_models)) {
        prev_model <- reuse_models[[m_idx]]
        is_jags_dsep    <- !is.null(prev_model$dsep_results)
        is_numpyro_dsep <- is.list(prev_model$dsep) && !is.null(prev_model$dsep$results)
        if (inherits(prev_model, "because") && (is_jags_dsep || is_numpyro_dsep)) {
            valid_models <- c(valid_models, list(prev_model))
        }
    }

    if (length(valid_models) == 0) {
        if (!quiet) {
            message("  No valid previous models found for reuse.")
        }
        return(list(found = list(), missing_indices = seq_along(current_tests)))
    }

    # Loop through current tests to find matches
    for (i in seq_along(current_tests)) {
        test_eq <- current_tests[[i]]

        # 1. Identify Required Support Equations for this specific test
        #    (Logic matches run_single_dsep_test extraction)
        req_eq_curr <- get_dsep_support_equations(
            test_eq,
            current_equations,
            family
        )

        # Normalize current test strings
        curr_test_str <- paste(deparse(test_eq), collapse = " ")
        curr_test_str <- gsub("\\s+", " ", trimws(curr_test_str))

        # Normalize required equations for comparison
        # We sort them to ensure order doesn't matter
        curr_req_strs <- sort(sapply(req_eq_curr, function(e) {
            str <- paste(deparse(e), collapse = " ")
            gsub("\\s+", " ", trimws(str))
        }))
        curr_req_hash <- paste(curr_req_strs, collapse = ";")

        found_match <- FALSE

        for (m_idx in seq_along(valid_models)) {
            prev_model <- valid_models[[m_idx]]

            # Prepare previous test strings
            prev_test_strs <- sapply(prev_model$dsep_tests, function(eq) {
                str <- paste(deparse(eq), collapse = " ")
                gsub("\\s+", " ", trimws(str))
            })

            match_idx <- match(curr_test_str, prev_test_strs)

            if (!is.na(match_idx)) {
                # Candidate found! Test formula matches.
                # CHECK: Do the required support equations match?
                prev_equations <- prev_model$equations %||%
                    prev_model$input$equations

                if (is.null(prev_equations)) {
                    # Cannot verify safety - skip
                    next
                }

                # Extract support equations from PREVIOUS model using PREVIOUS equations
                prev_family <- prev_model$input$family
                req_eq_prev <- get_dsep_support_equations(
                    test_eq,
                    prev_equations,
                    prev_family
                )

                prev_req_strs <- sort(sapply(req_eq_prev, function(e) {
                    str <- paste(deparse(e), collapse = " ")
                    gsub("\\s+", " ", trimws(str))
                }))
                prev_req_hash <- paste(prev_req_strs, collapse = ";")

                if (curr_req_hash == prev_req_hash) {
                    # MATCH CONFIRMED SAFE
                    # Retrieve from the correct slot depending on engine:
                    # JAGS/NIMBLE: $dsep_results is a list of coda-based result objects
                    # NumPyro:     $dsep$results is a data.frame (one row per test)
                    if (!is.null(prev_model$dsep_results)) {
                        found[[i]] <- prev_model$dsep_results[[match_idx]]
                    } else {
                        found[[i]] <- prev_model$dsep$results[match_idx, , drop = FALSE]
                    }
                    found_indices <- c(found_indices, i)
                    n_reused <- n_reused + 1
                    found_match <- TRUE
                    break # Stop searching models for this test; we found a reusable result
                }
            }
        }
    }

    missing_indices <- setdiff(seq_along(current_tests), found_indices)

    if (!quiet && n_reused > 0) {
        message(sprintf(
            "  Reusing %d d-separation tests from previous models.",
            n_reused
        ))
        if (length(missing_indices) == 0) {
            message("  All tests found in cache. Skipping new JAGS runs.")
        } else {
            message(sprintf(
                "  %d tests remain to be run.",
                length(missing_indices)
            ))
        }
    }

    return(list(found = found, missing_indices = missing_indices))
}

#' Helper to extract support equations (duplicates logic in run_single_dsep_test)
#' @keywords internal
get_dsep_support_equations <- function(test_eq, all_equations, family) {
    dsep_equations <- list(test_eq)

    if (!is.null(family) && any(family == "occupancy")) {
        added <- TRUE
        while (added) {
            added <- FALSE
            current_vars <- unique(unlist(lapply(dsep_equations, all.vars)))

            for (v in current_vars) {
                is_occ <- !is.na(family[v]) && family[v] == "occupancy"
                is_det <- grepl("^p_", v) &&
                    !is.na(family[sub("^p_", "", v)]) &&
                    family[sub("^p_", "", v)] == "occupancy"
                is_psi <- grepl("^psi_", v) &&
                    !is.na(family[sub("^psi_", "", v)]) &&
                    family[sub("^psi_", "", v)] == "occupancy"

                base_occ <- if (is_det) {
                    sub("^p_", "", v)
                } else if (is_psi) {
                    sub("^psi_", "", v)
                } else {
                    v
                }

                if (is_occ || is_det || is_psi) {
                    test_eq_resp <- as.character(test_eq)[2]

                    # Check for p_X and psi_X equations
                    targets <- c(
                        paste0("p_", base_occ),
                        paste0("psi_", base_occ)
                    )

                    for (tgt in targets) {
                        if (tgt != test_eq_resp) {
                            for (eq in all_equations) {
                                # Check if eq defines tgt
                                if (
                                    length(as.character(eq)) == 3 &&
                                        as.character(eq)[2] == tgt
                                ) {
                                    # Check if already in dsep_equations
                                    eq_str <- paste(deparse(eq), collapse = " ")
                                    exists <- any(sapply(
                                        dsep_equations,
                                        function(e) {
                                            paste(deparse(e), collapse = " ") ==
                                                eq_str
                                        }
                                    ))
                                    if (!exists) {
                                        dsep_equations <- c(
                                            dsep_equations,
                                            list(eq)
                                        )
                                        added <- TRUE
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return(dsep_equations)
}
