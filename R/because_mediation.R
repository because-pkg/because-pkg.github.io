#' Automated Mediation Analysis
#'
#' Decomposes the total effect of an exposure on an outcome into direct and indirect pathways
#' by tracing paths in the fitted DAG and multiplying posterior coefficients.
#'
#' @param fit A fitted object from \code{because()}.
#' @param exposure Character string; the name of the exposure variable.
#' @param outcome Character string; the name of the outcome variable.
#' @param prob Numeric; probability mass for the credible interval (default 0.95).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{paths}: A data frame summarizing each path (Path string, Mean, SD, CI).
#'   \item \code{summary}: A data frame summarizing Total, Direct, and Total Indirect effects.
#'   \item \code{samples}: A matrix of posterior samples for each path and the total effect.
#' }
#'
#' @details
#' This function reconstructs the causal graph from the \code{parameter_map} stored in the fit object.
#' It uses \code{igraph} to find all simple paths from \code{exposure} to \code{outcome}.
#' For each path, it identifies the corresponding regression coefficients (\code{beta_...})
#' and computes their product across all MCMC samples.
#'
#' \strong{Note:} This assumes linear relationships for indirect effects (\eqn{\beta_1 \times \beta_2}).
#' For non-linear models, this is an approximation of the average causal effect.
#'
#' @importFrom igraph graph_from_data_frame all_simple_paths as_ids
#' @importFrom stats sd quantile
#' @examples
#' \dontrun{
#' # Assuming a dataset 'df' exists with variables X, M, and Y
#' equations <- list(M ~ X, Y ~ M + X)
#' fit <- because(equations, data = df)
#' med_results <- because_mediation(fit, exposure = "X", outcome = "Y")
#' print(med_results$summary)
#' }
#'
#' @export
because_mediation <- function(fit, exposure, outcome, prob = 0.95) {
    if (is.null(fit$parameter_map)) {
        stop(
            "Fit object does not contain a parameter_map. Please refit the model with a newer version of 'because'."
        )
    }

    # 1. Build Graph from Parameter Map
    # The map contains: response, predictor, parameter
    # Since fit$parameter_map is a data frame, we can use it directly
    if (is.data.frame(fit$parameter_map)) {
        edges <- fit$parameter_map
        # Rename columns to match igraph expectations (from, to)
        colnames(edges)[colnames(edges) == "predictor"] <- "from"
        colnames(edges)[colnames(edges) == "response"] <- "to"
    } else {
        # Fallback if it is somehow a list (older versions?)
        edges <- do.call(
            rbind,
            lapply(fit$parameter_map, function(x) {
                data.frame(
                    from = x$predictor,
                    to = x$response,
                    parameter = x$parameter,
                    stringsAsFactors = FALSE
                )
            })
        )
    }

    # IMPORTANT: igraph uses the FIRST TWO columns as source/target.
    # We must ensure 'from' is 1st and 'to' is 2nd.
    if (!all(c("from", "to", "parameter") %in% colnames(edges))) {
        stop(
            "Parameter map structure unexpected. Columns 'from', 'to', 'parameter' required."
        )
    }

    # Reorder explicitly
    edges <- edges[, c("from", "to", "parameter")]

    # Filter to structural betas (exclude intercepts/auxiliary if any slipped in)
    # Ideally parameter_map only has structural predictors

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(
            "Package 'igraph' is required for mediation analysis. Please install it."
        )
    }

    g <- igraph::graph_from_data_frame(edges, directed = TRUE)

    # Check node existence
    if (!exposure %in% igraph::V(g)$name) {
        stop(paste("Exposure", exposure, "not found in model graph."))
    }
    if (!outcome %in% igraph::V(g)$name) {
        stop(paste("Outcome", outcome, "not found in model graph."))
    }

    # 2. Find All Simple Paths
    # Returns a list of vertex sequences
    paths_raw <- igraph::all_simple_paths(
        g,
        from = exposure,
        to = outcome,
        mode = "out"
    )

    if (length(paths_raw) == 0) {
        message("No directed paths found from ", exposure, " to ", outcome, ".")
        return(NULL)
    }

    # 3. Process Paths and Calculate Effects
    # Get MCMC samples
    samples_mat <- as.matrix(fit$samples) # Assuming fit$samples is the coda object
    if (is.null(samples_mat)) {
        # If samples not directly available, try to standard extraction
        stop("Could not extract MCMC samples from fit object.")
    }

    path_results <- list()
    total_effect_samples <- rep(0, nrow(samples_mat))
    direct_effect_samples <- rep(0, nrow(samples_mat))
    total_indirect_samples <- rep(0, nrow(samples_mat))

    path_names <- c()

    for (i in seq_along(paths_raw)) {
        p_nodes <- igraph::as_ids(paths_raw[[i]])
        # Construct path string: A -> B -> C
        p_str <- paste(p_nodes, collapse = " -> ")
        path_names[i] <- p_str

        # Identify edges and parameters along the path
        # Path is nodes: n1, n2, n3... Edge is (n1,n2), (n2,n3)
        path_effect <- rep(1, nrow(samples_mat))

        params_in_path <- c()

        for (j in 1:(length(p_nodes) - 1)) {
            u <- p_nodes[j]
            v <- p_nodes[j + 1]

            # Find parameter name in edges df
            # Should be unique for 'because' models (one edge per pair)
            param_name <- edges$parameter[edges$from == u & edges$to == v]

            if (length(param_name) == 0) {
                warning(paste(
                    "Parameter for edge",
                    u,
                    "->",
                    v,
                    "not found. Skipping path."
                ))
                path_effect <- 0
                break
            }

            # Check if sample exists
            # param_name might be an array beta[k] or something. Param map should have exact name.
            # If random effects formulation, structural betas are usually fixed effects.

            # Handle potential mismatched names in samples vs map?
            # Assuming map matches samples column names.
            if (!param_name %in% colnames(samples_mat)) {
                # Try to find matching column (e.g. if partial matching needed)
                # But 'because' usually keeps them clean.
                stop(paste(
                    "Parameter",
                    param_name,
                    "not found in MCMC samples."
                ))
            }

            params_in_path <- c(params_in_path, param_name)
            path_effect <- path_effect * samples_mat[, param_name]
        }

        # Store samples for this path
        path_results[[p_str]] <- path_effect

        # Accumulate Total Effect
        total_effect_samples <- total_effect_samples + path_effect

        # Classify as Direct or Indirect
        if (length(p_nodes) == 2) {
            # Direct path (Exposure -> Outcome)
            direct_effect_samples <- direct_effect_samples + path_effect
        } else {
            # Indirect path
            total_indirect_samples <- total_indirect_samples + path_effect
        }
    }

    # 4. Summarize Results
    summarize_vec <- function(vec, name) {
        m <- mean(vec)
        sd_val <- sd(vec)
        ci <- quantile(vec, probs = c((1 - prob) / 2, 1 - (1 - prob) / 2))
        data.frame(
            Type = name,
            Mean = m,
            SD = sd_val,
            Lower = unname(ci[1]),
            Upper = unname(ci[2]),
            stringsAsFactors = FALSE
        )
    }

    # Summary Table
    summary_df <- rbind(
        summarize_vec(total_effect_samples, "Total Effect"),
        summarize_vec(direct_effect_samples, "Direct Effect"),
        summarize_vec(total_indirect_samples, "Total Indirect Effect")
    )

    # Path Table
    path_df <- do.call(
        rbind,
        lapply(names(path_results), function(p) {
            res <- summarize_vec(path_results[[p]], "Path")
            transform(
                res,
                Path = p,
                Type = ifelse(
                    length(strsplit(p, " -> ")[[1]]) == 2,
                    "Direct",
                    "Indirect"
                )
            )
        })
    )
    # Clean up path df
    path_df <- path_df[, c("Path", "Type", "Mean", "SD", "Lower", "Upper")]

    # Compute Proportion Mediated (on likelihood of samples, not just means)
    # Prop = Indirect / Total
    # Caution: unstable if Total is near 0.
    prop_mediated_samples <- total_indirect_samples / total_effect_samples
    prop_summary <- summarize_vec(prop_mediated_samples, "Proportion Mediated")

    list(
        summary = summary_df,
        paths = path_df,
        proportion = prop_summary,
        samples = list(
            total = total_effect_samples,
            direct = direct_effect_samples,
            indirect = total_indirect_samples,
            paths = path_results
        )
    )
}
