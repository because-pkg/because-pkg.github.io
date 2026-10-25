#' Extract Imputed Values for Missing Data
#'
#' Extracts the posterior distributions of missing values that were imputed by the model.
#'
#' @param object A \code{because} model object.
#' @param id_col Optional character vector of IDs (e.g., species names) corresponding to the rows of the data.
#'   If the data in the model object already has names (e.g., from a phylogenetic model), these are used automatically.
#'   If provided, this vector must have the same length as the number of observations in the model (N).
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{Variable}: Name of the variable with missing data.
#'   \item \code{ID}: The identifier (e.g., Species) for the observation (if available).
#'   \item \code{RowIndex}: The original row index in the data.
#'   \item \code{Mean}: Posterior mean of the imputed value.
#'   \item \code{SD}: Posterior standard deviation.
#'   \item \code{Q2.5}: 2.5% quantile (lower credible interval).
#'   \item \code{Q50}: Median.
#'   \item \code{Q97.5}: 97.5% quantile (upper credible interval).
#' }
#'
#' @details
#' This function identifies which values in the original data were missing (\code{NA})
#' and looks up the corresponding imputed nodes in the posterior samples.
#'
#' Note: The imputed values are only available if they were monitored during the run.
#' You must use \code{monitor = "all"} (or include the variable names in \code{monitor})
#' when running \code{because()} to ensure these values are saved.
#'
#' @examples
#' \dontrun{
#' # Assuming 'fit' is a because model run with missing data and monitor="all"
#' imputed_values <- extract_imputed(fit)
#' head(imputed_values)
#' }
#'
#' @export
extract_imputed <- function(object, id_col = NULL) {
    if (!inherits(object, "because")) {
        stop("Object must be of class 'because'.")
    }

    data <- object$data
    # If data is not a list (e.g. legacy/unexpected), try to use as is provided it has names
    if (!is.list(data)) {
        stop("Model data structure not recognized (expected a list).")
    }

    # If id_col provided, validate length
    # We need N (sample size). JAGS data list usually has 'N'
    N <- data$N
    if (is.null(N)) {
        # Try to infer N from first vector variable
        vecs <- Filter(is.vector, data)
        if (length(vecs) > 0) N <- length(vecs[[1]])
    }

    # If id_col is not provided, try to use the one stored in the object
    if (is.null(id_col) && !is.null(object$species_order)) {
        id_col <- object$species_order
    }

    if (!is.null(id_col)) {
        if (!is.null(N) && length(id_col) != N) {
            stop(sprintf(
                "Provided id_col length (%d) does not match model sample size (%d).",
                length(id_col),
                N
            ))
        }
    }

    results_list <- list()

    # Iterate over variables in data
    for (var_name in names(data)) {
        values <- data[[var_name]]

        # Skip if not numeric/vector/matrix or if explicitly 'N' or 'ID'
        if (var_name %in% c("N", "ID", "zeros") || is.null(values)) {
            next
        }

        # Handle Vector (Variable[i])
        if (is.vector(values) && is.numeric(values)) {
            if (any(is.na(values))) {
                miss_idx <- which(is.na(values))

                for (i in miss_idx) {
                    # Construct parameter name: "Var[i]"
                    param_name <- paste0(var_name, "[", i, "]")

                    # Determine ID
                    row_id <- i
                    if (!is.null(id_col)) {
                        row_id <- id_col[i]
                    } else if (!is.null(names(values))) {
                        row_id <- names(values)[i]
                    }

                    stats <- get_stats(object, param_name)

                    if (!is.null(stats)) {
                        results_list[[length(results_list) + 1]] <- data.frame(
                            Variable = var_name,
                            ID = as.character(row_id),
                            RowIndex = i,
                            Mean = stats["Mean"],
                            SD = stats["SD"],
                            Q2.5 = stats["2.5%"],
                            Q50 = stats["50%"],
                            Q97.5 = stats["97.5%"],
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        } else if (is.matrix(values) && is.numeric(values)) {
            # Handle Matrix (Variable[i,j]) - e.g. reps
            if (any(is.na(values))) {
                miss_idx <- which(is.na(values), arr.ind = TRUE) # Returns matrix with row, col

                # Note: In 'because', usually matrix NAs in reps are just padded data or missing reps.
                # But 'because' models reps as 'Var_obs[i,j] ~ dnorm(...)'.
                # If the value is missing, it is imputed.
                # Parameter name: "Var[row,col]"

                for (k in seq_len(nrow(miss_idx))) {
                    r <- miss_idx[k, 1]
                    c <- miss_idx[k, 2]

                    param_name <- paste0(var_name, "[", r, ",", c, "]")

                    # Determine ID (row name)
                    row_id <- r
                    if (!is.null(id_col)) {
                        row_id <- id_col[r]
                    } else if (!is.null(rownames(values))) {
                        row_id <- rownames(values)[r]
                    }

                    stats <- get_stats(object, param_name)

                    if (!is.null(stats)) {
                        results_list[[length(results_list) + 1]] <- data.frame(
                            Variable = paste0(var_name, "_rep", c),
                            ID = as.character(row_id),
                            RowIndex = r,
                            Mean = stats["Mean"],
                            SD = stats["SD"],
                            Q2.5 = stats["2.5%"],
                            Q50 = stats["50%"],
                            Q97.5 = stats["97.5%"],
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
    }

    if (length(results_list) == 0) {
        message(
            "No imputed values found. Either no data was missing, or imputed nodes were not monitored (did you use monitor='all'?)."
        )
        return(data.frame())
    }

    do.call(rbind, results_list)
}

# Helper to extract stats from summary or samples
get_stats <- function(object, param_name) {
    # Try summary first
    if (!is.null(object$summary)) {
        sum_stats <- object$summary
        if (param_name %in% rownames(sum_stats$statistics)) {
            mean_val <- sum_stats$statistics[param_name, "Mean"]
            sd_val <- sum_stats$statistics[param_name, "SD"]
            q2.5 <- sum_stats$quantiles[param_name, "2.5%"]
            q50 <- sum_stats$quantiles[param_name, "50%"]
            q97.5 <- sum_stats$quantiles[param_name, "97.5%"]
            return(c(
                "Mean" = mean_val,
                "SD" = sd_val,
                "2.5%" = q2.5,
                "50%" = q50,
                "97.5%" = q97.5
            ))
        }
    }

    # If not in summary, try samples (calculate on fly) - slower but robust
    # Check if param_name is in samples columns
    if (!is.null(object$samples)) {
        # Combine chains
        # Usually samples is an mcmc.list
        varnames <- coda::varnames(object$samples)
        if (param_name %in% varnames) {
            # Extract data
            mcmc_mat <- as.matrix(object$samples[, param_name])
            mean_val <- mean(mcmc_mat)
            sd_val <- sd(mcmc_mat)
            quants <- quantile(mcmc_mat, probs = c(0.025, 0.5, 0.975))
            return(c(
                "Mean" = mean_val,
                "SD" = sd_val,
                "2.5%" = quants[1],
                "50%" = quants[2],
                "97.5%" = quants[3]
            ))
        }
    }

    return(NULL)
}
