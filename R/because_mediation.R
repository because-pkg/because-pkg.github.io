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
    if (!inherits(fit, "because") || is.null(fit$parameter_map)) {
        stop(
            "'fit' must be a valid 'because' model object containing a 'parameter_map'."
        )
    }

    # 1. Build Graph from Parameter Map
    edges <- fit$parameter_map
    colnames(edges)[colnames(edges) == "predictor"] <- "from"
    colnames(edges)[colnames(edges) == "response"] <- "to"

    if (!all(c("from", "to", "parameter") %in% colnames(edges))) {
        stop(
            "Parameter map structure unexpected. Columns 'from', 'to', 'parameter' required."
        )
    }

    # igraph::graph_from_data_frame treats the first two columns as source -> target
    edges <- edges[, c("from", "to", "parameter")]

    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop(
            "Package 'igraph' is required for mediation analysis. Please install it."
        )
    }

    g <- igraph::graph_from_data_frame(edges, directed = TRUE)

    if (!exposure %in% igraph::V(g)$name) {
        stop(paste("Exposure", exposure, "not found in model graph."))
    }
    if (!outcome %in% igraph::V(g)$name) {
        stop(paste("Outcome", outcome, "not found in model graph."))
    }

    # 2. Find All Simple Paths
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
    samples_mat <- as.matrix(fit$samples)
    if (is.null(samples_mat)) {
        stop("Could not extract MCMC samples from fit object.")
    }

    path_results <- list()
    total_effect_samples <- rep(0, nrow(samples_mat))
    direct_effect_samples <- rep(0, nrow(samples_mat))
    total_indirect_samples <- rep(0, nrow(samples_mat))

    path_names <- c()

    for (i in seq_along(paths_raw)) {
        p_nodes <- igraph::as_ids(paths_raw[[i]])
        p_str <- paste(p_nodes, collapse = " -> ")
        path_names[i] <- p_str

        path_effect <- rep(1, nrow(samples_mat))
        params_in_path <- c()

        for (j in 1:(length(p_nodes) - 1)) {
            u <- p_nodes[j]
            v <- p_nodes[j + 1]

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

            if (!param_name %in% colnames(samples_mat)) {
                stop(paste(
                    "Parameter",
                    param_name,
                    "not found in MCMC samples."
                ))
            }

            params_in_path <- c(params_in_path, param_name)
            path_effect <- path_effect * samples_mat[, param_name]
        }

        path_results[[p_str]] <- path_effect
        total_effect_samples <- total_effect_samples + path_effect

        if (length(p_nodes) == 2) {
            direct_effect_samples <- direct_effect_samples + path_effect
        } else {
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

    summary_df <- rbind(
        summarize_vec(total_effect_samples, "Total Effect"),
        summarize_vec(direct_effect_samples, "Direct Effect"),
        summarize_vec(total_indirect_samples, "Total Indirect Effect")
    )

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
    path_df <- path_df[, c("Path", "Type", "Mean", "SD", "Lower", "Upper")]

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
