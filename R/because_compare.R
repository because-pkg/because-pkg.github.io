#' Compare Because Models
#'
#' A unified function to either (1) compare previously fitted models, or (2) run multiple model specifications in parallel and then compare them.
#'
#' @param ... For comparing fitted models: individual fitted model objects of class \code{"because"}.
#'   For running models: additional arguments passed to \code{\link{because}} (e.g., \code{n.iter}).
#' @param model_specs A named list of model specifications to run (Mode 2). Each element should be a list containing arguments for \code{because}.
#'   Alternatively, this argument can accept the first fitted model object (Mode 1).
#' @param data The dataset (required for Mode 2). Alternatively, the second fitted model object (Mode 1).
#' @param tree The phylogenetic tree (optional for Mode 2). Alternatively, the third fitted model object (Mode 1).
#' @param n.cores Number of cores for parallel execution (Mode 2). Default is 1.
#' @param cl Optional cluster object (Mode 2).
#' @param sort Logical. If \code{TRUE} (default), sort comparison table by WAIC.
#'
#' @return
#' If comparing fitted models: A class \code{"because_comparison"} object (data frame) with WAIC rankings.
#'
#' If running models: A list containing:
#' \item{results}{List of fitted model objects.}
#' \item{comparison}{The comparison data frame.}
#'
#' @details
#' \strong{Mode 1: Compare Fitted Models}
#' Call \code{because_compare(fit1, fit2)} or \code{because_compare(models = list(fit1, fit2))}.
#' Extracts WAIC (with SE) from each model and ranks them.
#'
#' \strong{Mode 2: Run and Compare}
#' Call \code{because_compare(model_specs = list(m1=..., m2=...), data=data, tree=tree)}.
#' This runs the models in parallel and returns the comparison.
#'
#' @examples
#' \dontrun{
#'   # Mode 1: Compare existing fits
#'   because_compare(fit1, fit2)
#'
#'   # Mode 2: Run and compare
#'   specs <- list(m1 = list(equations = list(Y ~ X)), m2 = list(equations = list(Y ~ X + Z)))
#'   res <- because_compare(specs, data = df, tree = tr, n.cores = 2)
#'   print(res$comparison)
#' }
#'
#' @export
because_compare <- function(
    ...,
    model_specs = NULL,
    data = NULL,
    tree = NULL,
    n.cores = 1,
    cl = NULL,
    sort = TRUE
) {
    # --- ARGUMENT HARVESTING & MODE DETECTION ---

    # Collect all named arguments and ... arguments
    dots <- list(...)

    # Check if we are in Mode 2 (Run Specs)
    # Criteria: 'model_specs' is a list of specs, OR dots[[1]] is a list of specs

    specs <- NULL
    run_data <- NULL
    run_tree <- NULL

    is_spec_list <- function(x) {
        is.list(x) &&
            !is.null(names(x)) &&
            all(sapply(x, function(s) {
                is.list(s) &&
                    ("equations" %in% names(s) || "formula" %in% names(s))
            }))
    }

    # Check explicit argument
    if (is_spec_list(model_specs)) {
        specs <- model_specs
        run_data <- data
        run_tree <- tree
    } else if (length(dots) >= 1 && is_spec_list(dots[[1]])) {
        # Check first positional argument (if model_specs was NULL)
        specs <- dots[[1]]

        # Only pull data/tree from UNNAMED positional arguments
        # names(dots) for unnamed is ""
        unnamed_indices <- which(names(dots) == "" | is.null(names(dots)))
        # Filter out specs itself (index 1)
        unnamed_indices <- setdiff(unnamed_indices, 1)

        # Attempt to retrieve data/tree from dots or named args
        if (!is.null(data)) {
            run_data <- data
        } else if (length(unnamed_indices) >= 1) {
            run_data <- dots[[unnamed_indices[1]]]
        }

        if (!is.null(tree)) {
            run_tree <- tree
        } else if (length(unnamed_indices) >= 2) {
            run_tree <- dots[[unnamed_indices[2]]]
        }
    }

    # If specs found, Execute Mode 2
    if (!is.null(specs)) {
        if (is.null(run_data)) {
            stop("Argument 'data' is required for running models.")
        }

        # Extract extra args for because from dots (excluding rescued pos args)
        # Handling mixed positional and named arguments in ...
        common_args <- dots[nzchar(names(dots))]

        # Run Function
        run_model_internal <- function(name, spec, d, t, extra) {
            # Remove WAIC and quiet from extra to avoid duplicates (hardcoded below)
            extra$WAIC <- NULL
            extra$quiet <- NULL

            args <- c(
                list(data = d, structure = t, WAIC = TRUE, quiet = TRUE),
                spec,
                extra
            )
            tryCatch(
                {
                    do.call(because::because, args)
                },
                error = function(e) {
                    message(paste("Model", name, "failed:", e$message))
                    warning(paste("Model", name, "failed:", e$message))
                    return(NULL)
                }
            )
        }

        message(sprintf(
            "Running %d models in comparison set...",
            length(specs)
        ))

        fit_results <- list()
        spec_names <- names(specs)

        # Handle if n.cores was passed as a model object by mistake
        # Ensure n.cores is numeric
        n_cores_val <- if (is.numeric(n.cores)) n.cores else 1

        if (n_cores_val > 1) {
            created_cluster <- FALSE
            if (is.null(cl)) {
                cl <- parallel::makeCluster(n_cores_val)
                on.exit(parallel::stopCluster(cl), add = TRUE)
                created_cluster <- TRUE
            }
            parallel::clusterEvalQ(cl, library(because))

            fit_results <- parallel::parLapply(cl, spec_names, function(nm) {
                run_model_internal(
                    nm,
                    specs[[nm]],
                    run_data,
                    run_tree,
                    common_args
                )
            })
            names(fit_results) <- spec_names
        } else {
            fit_results <- lapply(spec_names, function(nm) {
                message(paste("Fitting model:", nm))
                run_model_internal(
                    nm,
                    specs[[nm]],
                    run_data,
                    run_tree,
                    common_args
                )
            })
            names(fit_results) <- spec_names
        }

        # Filter failed
        failed <- sapply(fit_results, is.null)
        if (any(failed)) {
            message("Some models failed.")
        }
        fit_results <- fit_results[!failed]

        if (length(fit_results) == 0) {
            return(NULL)
        }

        # Compare
        comp <- because_compare(models = fit_results, sort = sort)

        return(list(results = fit_results, comparison = comp))
    } else {
        # --- Mode 1: Compare Fitted Models ---

        # Harvest models from ALL arguments
        candidates <- list(model_specs, data, tree, n.cores, cl) # Potential positional args
        candidates <- c(candidates, dots)

        # If user passed a list of models explicitly
        if (!is.null(dots$models) && is.list(dots$models)) {
            candidates <- c(candidates, dots$models)
        }

        # Filter for because objects
        is_fit <- function(x) inherits(x, "because")
        models <- Filter(is_fit, candidates)

        # Name models if needed
        if (length(models) > 0) {
            # Separate named vs unnamed arguments
            # Use names if available in source list, else generate
            curr_names <- names(models)
            if (is.null(curr_names)) {
                curr_names <- rep("", length(models))
            }

            for (i in seq_along(models)) {
                if (curr_names[i] == "") curr_names[i] <- paste0("Model_", i)
            }
            names(models) <- curr_names
        } else {
            stop(
                "No fitted 'because' models provided (Mode 1), nor valid 'model_specs' (Mode 2)."
            )
        }

        return(calc_waic_table(models, sort))
    }
}

# Internal Helper
calc_waic_table <- function(models, sort = TRUE) {
    res_list <- list()
    pointwise_list <- list()

    for (name in names(models)) {
        mod <- models[[name]]
        if (is.null(mod$WAIC)) {
            warning(paste("Model", name, "missing WAIC. Skipping."))
            next
        }

        res_list[[name]] <- c(
            WAIC = mod$WAIC["waic", "Estimate"],
            SE = mod$WAIC["waic", "SE"],
            p_waic = mod$WAIC["p_waic", "Estimate"]
        )

        pt <- attr(mod$WAIC, "pointwise")
        if (!is.null(pt)) pointwise_list[[name]] <- pt
    }

    if (length(res_list) == 0) {
        return(NULL)
    }

    df <- as.data.frame(do.call(rbind, res_list))
    df$Model <- rownames(df)

    if (sort) {
        df <- df[order(df$WAIC), ]
        pointwise_list <- pointwise_list[df$Model]
    }

    df$dWAIC <- df$WAIC - df$WAIC[1]

    dSE <- numeric(nrow(df))
    best_pt <- pointwise_list[[1]]

    for (i in seq_len(nrow(df))) {
        if (i == 1) {
            dSE[i] <- 0
        } else {
            curr_pt <- pointwise_list[[i]]
            if (
                !is.null(best_pt) &&
                    !is.null(curr_pt) &&
                    nrow(best_pt) == nrow(curr_pt) &&
                    "waic" %in% colnames(best_pt)
            ) {
                diffs <- curr_pt[, "waic"] - best_pt[, "waic"]
                dSE[i] <- sqrt(length(diffs) * var(diffs))
            } else {
                dSE[i] <- NA
            }
        }
    }
    df$dSE <- dSE

    rel_lik <- exp(-0.5 * df$dWAIC)
    df$weight <- rel_lik / sum(rel_lik)

    df <- df[, c("WAIC", "SE", "dWAIC", "dSE", "p_waic", "weight")]
    class(df) <- c("because_comparison", "data.frame")
    return(df)
}

#' @keywords internal
#' @export
print.because_comparison <- function(x, digits = 2, ...) {
    cat("Model Comparison (ordered by WAIC):\n")
    print.data.frame(x, digits = digits, ...)
    cat("\n")
}
