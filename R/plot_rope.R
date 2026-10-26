#' Plot Region of Practical Equivalence (ROPE)
#'
#' Generic for plotting posterior density distributions and visualizing their overlap 
#' with a Region of Practical Equivalence (ROPE).
#'
#' @param object A fitted model object.
#' @param ... Additional arguments.
#' @export
plot_rope <- function(object, ...) {
    UseMethod("plot_rope")
}

utils::globalVariables(c("Value", "x"))

#' @rdname plot_rope
#'
#' @param object A `because` object.
#' @param rope Optional ROPE specification. Can be \code{NULL} or \code{"default"} (uses principled
#'   family defaults: \code{c(-0.10, 0.10)} for Gaussian, \code{c(-0.05, 0.05)} for Poisson/NegBin,
#'   \code{c(-0.18, 0.18)} for Binomial/Ordinal/Multinomial), a numeric vector of length 2
#'   (e.g., \code{c(-0.1, 0.1)} applied globally), or a named list specifying custom intervals by family
#'   (e.g., \code{list(poisson = c(-0.05, 0.05))}), by response variable, or by parameter name.
#' @param parameters Optional character vector to filter which parameters to plot, or `"all"` 
#'   to plot all monitored parameters. If `NULL` (default), only structural path coefficients (`beta_*`) 
#'   and correlations (`rho_*`) are plotted for standard models, excluding nuisance parameters (intercepts 
#'   and variances).
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' The Region of Practical Equivalence (ROPE) formalizes Bayesian equivalence testing
#' (Kruschke 2018; Makowski et al. 2019). When evaluating causal paths or conditional
#' independence claims in d-separation, failure to reject a sharp null (\eqn{\beta = 0})
#' does not provide evidence of absence. ROPE defines an explicit negligible boundary
#' \eqn{[-\delta, +\delta]} around zero:
#' \itemize{
#'   \item \strong{Gaussian}: The default interval is \code{c(-0.10, 0.10)}, representing
#'     half of Cohen's negligible effect size benchmark (\eqn{\pm 0.1 \times \text{SD}})
#'     for standardized continuous traits (Kruschke 2018).
#'   \item \strong{Poisson & Negative Binomial}: Coefficients sit on the log-rate scale
#'     (\eqn{\log \text{IRR}}). The default interval \code{c(-0.05, 0.05)} represents an
#'     approximate \eqn{\pm 5\%} relative change in expected count
#'     (\eqn{\text{IRR} \in [\exp(-0.05), \exp(0.05)] \approx [0.95, 1.05]}).
#'   \item \strong{Binomial, Ordinal, & Multinomial}: In logistic models, the binary outcome
#'     is modeled via an underlying continuous latent variable with standard logistic error variance
#'     \eqn{\sigma^2 = \pi^2 / 3} and standard deviation \eqn{\sigma = \pi / \sqrt{3} \approx 1.814}.
#'     Applying Kruschke's \eqn{0.1 \times \sigma} rule to this latent scale yields
#'     \eqn{\pm 0.1 \times (\pi / \sqrt{3}) \approx \pm 0.181} (Gelman et al. 2008; Makowski et al. 2019).
#'     This corresponds to an Odds Ratio window of \eqn{\text{OR} \in [\exp(-0.181), \exp(0.181)] \approx [0.83, 1.20]}.
#' }
#'
#' @references
#' Kruschke, J. K. (2018). Rejecting or accepting parameter values in Bayesian estimation.
#' \emph{Advances in Methods and Practices in Psychological Science}, 1(2), 270–280. \doi{10.1177/2515245918771304}
#'
#' Makowski, D., Ben-Shachar, M. S., & Lüdecke, D. (2019). bayestestR: Describing effects
#' and their uncertainty, existence and significance within the Bayesian framework.
#' \emph{Journal of Open Source Software}, 4(40), 1541. \doi{10.21105/joss.01541}
#'
#' Gelman, A., Jakulin, A., Pittau, M. G., & Su, Y. S. (2008). A weakly informative
#' default prior distribution for logistic and other regression models.
#' \emph{The Annals of Applied Statistics}, 2(4), 1360–1383. \doi{10.1214/08-AOAS191}
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_manual geom_vline geom_segment labs theme_minimal theme element_text after_stat
plot_rope.because <- function(object, rope = "default", parameters = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required. Please install it.")
    }
    if (!requireNamespace("ggridges", quietly = TRUE)) {
        stop("Package 'ggridges' is required for plot_rope. Please install it.")
    }

    is_dsep <- !is.null(object$dsep) && (isTRUE(object$dsep) || is.list(object$dsep))
    
    plot_data <- data.frame()
    
    if (is_dsep) {
        # Check for pre-computed numpyro dsep results
        if (is.list(object$dsep) && !is.null(object$dsep$results)) {
            stop("ROPE plots are not currently supported for numpyro d-sep tests because the raw samples are not preserved.")
        }
        
        summ <- summary(object)
        res_table <- summ$results
        
        if (is.null(res_table) || nrow(res_table) == 0) {
            stop("No d-separation test results found.")
        }
        
        dsep_results <- object$dsep_results
        
        for (i in seq_len(nrow(res_table))) {
            test_label <- res_table$Test[i]
            param_name <- res_table$Parameter[i]
            
            if (!is.null(parameters) && !(param_name %in% parameters || test_label %in% parameters)) {
                next
            }
            
            # Find the samples in dsep_results
            for (dr in dsep_results) {
                if (!is.null(dr)) {
                    samps <- as.matrix(dr$samples)
                    if (param_name %in% colnames(samps)) {
                        # Simplify test label by removing | {}
                        clean_label <- gsub(" \\| \\{.*?\\}", "", test_label)
                        clean_label <- trimws(clean_label)
                        
                        # Infer response and family
                        resp_name <- sub("\\s*_\\|\\|_.*", "", test_label)
                        fam_i <- if (!is.null(object$family) && resp_name %in% names(object$family)) object$family[[resp_name]] else "gaussian"
                        r_int <- resolve_rope_for_param(param_name, resp = resp_name, fam = fam_i, rope_spec = rope)

                        vals <- samps[, param_name]
                        tmp <- data.frame(
                            Label = clean_label,
                            Value = vals,
                            Response = resp_name,
                            Family = fam_i,
                            r_min = r_int[1],
                            r_max = r_int[2],
                            stringsAsFactors = FALSE
                        )
                        plot_data <- rbind(plot_data, tmp)
                        break
                    }
                }
            }
        }
    } else {
        # Standard model
        summ <- summary(object)
        target_params <- rownames(summ$results)
        
        if (!is.null(parameters)) {
            if (!identical(parameters, "all")) {
                target_params <- intersect(target_params, parameters)
            }
        } else {
            # By default, focus on structural paths (betas) and correlations (rhos),
            # matching plot_coef(), excluding intercepts (alpha) and variances (sigma/tau)
            target_params <- target_params[grepl("^(beta_|rho_)", target_params)]
            target_params <- target_params[!grepl("(_det_|_1_times_)", target_params)]
        }
        
        if (length(target_params) == 0) {
            stop("No valid parameters found to plot.")
        }
        
        all_samps <- as.matrix(object$samples)
        
        for (p in target_params) {
            if (p %in% colnames(all_samps)) {
                vals <- all_samps[, p]

                # Resolve friendly label, response, and family
                label <- p
                resp <- "All"
                fam_i <- "gaussian"
                if (!is.null(object$parameter_map)) {
                    p_base <- gsub("\\[\\d+\\]$", "", p)
                    p_match <- object$parameter_map[object$parameter_map$parameter == p_base |
                                                   paste0(object$parameter_map$parameter, "[]") == p_base, ]
                    if (nrow(p_match) > 0) {
                        resp <- p_match$response[1]
                        pred <- p_match$predictor[1]
                        idx_match <- regmatches(p, regexpr("\\[(\\d+)\\]", p))
                        if (length(idx_match) > 0) {
                            idx <- gsub("[\\[\\]]", "", idx_match)
                            label <- paste0(resp, " ~ ", pred, " (", idx, ")")
                        } else {
                            label <- paste0(resp, " ~ ", pred)
                        }
                    }
                }
                if (!is.null(resp) && resp != "All" && !is.null(object$family) && resp %in% names(object$family)) {
                    fam_i <- object$family[[resp]]
                }

                r_int <- resolve_rope_for_param(p, resp = resp, fam = fam_i, rope_spec = rope)

                tmp <- data.frame(
                    Label = label,
                    Value = vals,
                    Response = resp,
                    Family = fam_i,
                    r_min = r_int[1],
                    r_max = r_int[2],
                    stringsAsFactors = FALSE
                )
                plot_data <- rbind(plot_data, tmp)
            }
        }
    }
    
    if (nrow(plot_data) == 0) {
        stop("No data available to plot.")
    }
    
    # Ensure Label is a factor to preserve order (bottom to top in ggridges)
    # Reverse so the first test is at the top
    plot_data$Label <- factor(plot_data$Label, levels = rev(unique(plot_data$Label)))
    
    unique_labels <- levels(plot_data$Label)
    r_mins_by_level <- numeric(length(unique_labels))
    r_maxs_by_level <- numeric(length(unique_labels))
    for (lvl_idx in seq_along(unique_labels)) {
        lbl <- unique_labels[lvl_idx]
        sub_row <- plot_data[plot_data$Label == lbl, ]
        r_mins_by_level[lvl_idx] <- sub_row$r_min[1]
        r_maxs_by_level[lvl_idx] <- sub_row$r_max[1]
    }
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, y = Label)) +
        ggridges::stat_density_ridges(
            geom = "density_ridges_gradient",
            ggplot2::aes(fill = ggplot2::after_stat(x >= r_mins_by_level[group] & x <= r_maxs_by_level[group])),
            color = "black",
            scale = 0.95,
            rel_min_height = 0.005,
            show.legend = TRUE
        ) +
        ggplot2::scale_fill_manual(
            name = "Region",
            values = c("TRUE" = "#89C5C5", "FALSE" = "#F8766D"),
            labels = c(
                "TRUE" = "Inside ROPE (Equivalence)", 
                "FALSE" = "Outside ROPE"
            ),
            breaks = c("TRUE", "FALSE")
        )
        
    all_same_rope <- (length(unique(r_mins_by_level)) == 1) && (length(unique(r_maxs_by_level)) == 1)
    if (all_same_rope) {
        p <- p + ggplot2::geom_vline(xintercept = c(r_mins_by_level[1], r_maxs_by_level[1]), linetype = "dashed", color = "darkgray", linewidth = 1)
    } else {
        # Draw per-ridge vertical boundary segments
        v_lines_df <- data.frame(
            x = c(r_mins_by_level, r_maxs_by_level),
            y = rep(seq_along(unique_labels), 2),
            yend = rep(seq_along(unique_labels) + 0.85, 2)
        )
        # Draw horizontal baseline bracket indicating the ROPE interval
        h_lines_df <- data.frame(
            x = r_mins_by_level,
            xend = r_maxs_by_level,
            y = seq_along(unique_labels)
        )
        p <- p +
            ggplot2::geom_segment(
                data = v_lines_df,
                ggplot2::aes(x = x, xend = x, y = y, yend = yend),
                linetype = "dashed",
                color = "darkgray",
                linewidth = 0.8,
                inherit.aes = FALSE
            ) +
            ggplot2::geom_segment(
                data = h_lines_df,
                ggplot2::aes(x = x, xend = xend, y = y, yend = y),
                color = "#4A8F8F",
                linewidth = 1.2,
                inherit.aes = FALSE
            )
    }

    # Format human-readable explanation of active ROPE specifications
    rope_info <- explain_rope_spec(unique(plot_data[, c("Response", "Family", "r_min", "r_max")]), data = object$data %||% object$original_data)
    sub_title_text <- if (length(rope_info$lines) > 0) paste(rope_info$lines, collapse = " | ") else "Density curves stacked vertically"
    if (!is.null(rope_info$note)) {
        sub_title_text <- paste0(sub_title_text, "\n", rope_info$note)
    }

    p <- p +
        ggplot2::labs(
            title = if(is_dsep) "d-separation Tests: ROPE Overlap" else "Posterior Distributions: ROPE Overlap",
            subtitle = sub_title_text,
            x = "Effect Size",
            y = ""
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = "bottom"
        )
        
    return(p)
}
