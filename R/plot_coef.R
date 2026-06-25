#' @rdname plot_coef
#'
#' @description 
#' Creates a 'caterpillar plot' (point and whisker) of all path coefficients 
#' from a `because` model. This visualization provides a precise statistical 
#' complement to the `plot_dag()` overview, allowing for a side-by-side 
#' comparison of effect sizes and their Bayesian credibility intervals.
#'
#' @param object A `because` object.
#' @param type Character; either `"raw"` (default) or `"marginal"`.
#'   - `"marginal"`: Shows **Average Marginal Effects (AME)**. For categorical 
#'     predictors, this represents the average shift in the outcome (e.g. probability or counts) 
#'     associated with a one-category change. This is the recommended scale for 
#'     comparing cross-model impacts.
#'   - `"raw"` (default): Shows the raw structural parameters (betas/rhos) from the JAGS model. 
#'     Useful for model diagnostics but harder to interpret on the original data scale.
#' @param multinomial_probabilities Logical; if `TRUE` (default), expands 
#'   multinomial predictors into a "bundle" of category-specific effects, 
#'   matching the arcs in the DAG.
#' @param color_scheme Character; color scheme for significance. Options:
#'   - `"sig_only"` (default): Discrete Black/Grey scheme. Significant paths 
#'     (where 95% CI excludes zero) are Black; non-significant are Light Grey.
#'   - `"directional"`: Switched to a directional Red/Blue/Grey scheme.
#'   - `"monochrome"`: All effects are rendered in Black regardless of significance.
#' @param prob Numeric; probability mass for the credibility interval (default 0.95).
#' @param ... Additional arguments.
#'
#' @return A `ggplot` object. Use standard `ggplot2` functions like `+ ggtitle()` 
#'   or `+ theme()` to further customize the output.
#'
#' @details 
#' The function automatically sorts coefficients by the **Response Variable**, 
#' effectively grouping all predictors for a given outcome together. This 
#' hierarchy makes it intuitive to read 'down' the plot to see what factors 
#' contribute most to a specific part of your causal system.
#'
#' @examples
#' \dontrun{
#' # Fit a model
#' fit <- because(list(Y ~ X + Z, X ~ Z), family = c(Y="binomial", X="gaussian"), data = dat)
#' 
#' # Plot marginal effects for intuitive interpretation
#' plot_coef(fit, type = "marginal")
#' 
#' # Customizing the plot
#' library(ggplot2)
#' plot_coef(fit) + labs(title = "Causal Influence on Wildlife Tolerance")
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_hline coord_flip labs theme_minimal scale_color_manual
plot_coef.because <- function(
    object,
    type = "raw",
    multinomial_probabilities = TRUE,
    color_scheme = "sig_only",
    prob = 0.95,
    ...
) {
    if (!inherits(object, "because")) {
        stop("object must be of class 'because'")
    }

    # 1. Extract Data
    me_table <- NULL
    stats <- NULL
    quantiles <- NULL

    if (type == "marginal") {
        if (is.null(object$parameter_map)) {
            stop("Marginal effects require a parameter map. Did you fit the model with because()?")
        }
        # Use high precision for all marginal plots
        set.seed(12345)
        me_table <- marginal_effects(
            object, 
            samples = 1000, 
            prob = prob,
            multinomial_probabilities = multinomial_probabilities
        )
    } else {
        if (prob == 0.95 && !is.null(object$summary)) {
            stats <- object$summary$statistics
            quantiles <- object$summary$quantiles
        } else {
            lower_prob <- (1 - prob) / 2
            upper_prob <- 1 - lower_prob
            quantiles_req <- c(lower_prob, 0.5, upper_prob)
            summ <- summary(object$samples, quantiles = quantiles_req)
            stats <- summ$statistics
            quantiles <- summ$quantiles
        }
    }

    # 2. Build Plotting Data Frame
    plot_df <- data.frame(
        Path = character(),
        Estimate = numeric(),
        Lower = numeric(),
        Upper = numeric(),
        Significant = character(),
        Response = character(),
        stringsAsFactors = FALSE
    )

    if (type == "marginal") {
        # Process Marginal Effects Table
        for (i in seq_len(nrow(me_table))) {
            resp <- as.character(me_table$Response[i])
            pred <- as.character(me_table$Predictor[i])
            cat  <- me_table$Category[i]
            est  <- me_table$Effect[i]
            low  <- me_table$Lower[i]
            upp  <- me_table$Upper[i]

            path_label <- if (!is.na(cat)) {
                paste0(resp, " ~ ", pred, " (", cat, ")")
            } else {
                paste0(resp, " ~ ", pred)
            }

            sig <- "non-sig"
            if (sign(low) == sign(upp)) {
                if (color_scheme == "directional") {
                    sig <- if (est > 0) "pos" else "neg"
                } else {
                    sig <- "sig"
                }
            }

            plot_df <- rbind(plot_df, data.frame(
                Path = path_label,
                Estimate = est,
                Lower = low,
                Upper = upp,
                Significant = sig,
                Response = resp,
                stringsAsFactors = FALSE
            ))
        }
    } else {
         # Process Raw Parameters
         # We focus on betas (paths) and rhos (correlations)
         param_names <- rownames(stats)
         target_params <- param_names[grepl("^(beta_|rho_)", param_names)]
         
         # Exclude internal deterministic link parameters
         target_params <- target_params[!grepl("(_det_|_1_times_)", target_params)]

         # Helpers for dynamically resolving the quantile column names
         format_pct <- function(p) paste0(format(100 * p, trim = TRUE, digits = max(2, getOption("digits"))), "%")
         lower_prob <- (1 - prob) / 2
         upper_prob <- 1 - lower_prob
         lower_col <- format_pct(lower_prob)
         upper_col <- format_pct(upper_prob)

         for (pname in target_params) {
             est  <- stats[pname, "Mean"]
             
             # Safely get lower bound
             if (lower_col %in% colnames(quantiles)) {
                 low  <- quantiles[pname, lower_col]
             } else {
                 low  <- quantiles[pname, 1]
             }
             
             # Safely get upper bound
             if (upper_col %in% colnames(quantiles)) {
                 upp  <- quantiles[pname, upper_col]
             } else {
                 upp  <- quantiles[pname, ncol(quantiles)]
             }

             # Try to resolve friendly names from parameter_map
             label <- pname
             resp_grp <- "All"
             if (!is.null(object$parameter_map)) {
                 # Handle arrays (e.g. beta_Y_X[2]) vs scalars (beta_Y_X)
                 p_base <- gsub("\\[\\d+\\]$", "", pname)
                 p_match <- object$parameter_map[object$parameter_map$parameter == p_base | 
                                               paste0(object$parameter_map$parameter, "[]") == p_base, ]
                 if (nrow(p_match) > 0) {
                     resp <- p_match$response[1]
                     pred <- p_match$predictor[1]
                     resp_grp <- resp
                     
                     # Check if it's an array index (category)
                     idx_match <- regmatches(pname, regexpr("\\[(\\d+)\\]", pname))
                     if (length(idx_match) > 0) {
                         idx <- gsub("[\\[\\]]", "", idx_match)
                         label <- paste0(resp, " ~ ", pred, " (", idx, ")")
                     } else {
                         label <- paste0(resp, " ~ ", pred)
                     }
                 }
             }

             sig <- "non-sig"
             if (sign(low) == sign(upp)) {
                 if (color_scheme == "directional") {
                     sig <- if (est > 0) "pos" else "neg"
                 } else {
                     sig <- "sig"
                 }
             }

             plot_df <- rbind(plot_df, data.frame(
                 Path = label,
                 Estimate = est,
                 Lower = low,
                 Upper = upp,
                 Significant = sig,
                 Response = resp_grp,
                 stringsAsFactors = FALSE
             ))
         }
    }

    if (nrow(plot_df) == 0) {
        stop("No causal paths or correlations found to plot.")
    }

    # 3. Rendering
    # Sort: Alphabetical by Response (group), and then alphabetical by Path
    plot_df <- plot_df[order(plot_df$Response, plot_df$Path), ]
    plot_df$Path <- factor(plot_df$Path, levels = rev(unique(plot_df$Path)))

    color_map <- c(
        "pos" = "#000000",     # Black
        "neg" = "#000000",     # Black
        "sig" = "#000000",     # Black
        "non-sig" = "#bdbdbd"  # Light Grey
    )

    # Symmetric axis: always show both sides of zero so positive and negative
    # effects can be visually compared regardless of the data range.
    # Use coord_flip(ylim=) rather than scale_y_continuous(limits=) to zoom
    # the coordinate system without clipping whiskers at the scale level.
    max_abs <- max(abs(c(plot_df$Lower, plot_df$Upper, plot_df$Estimate)), na.rm = TRUE)
    axis_lim <- c(-max_abs, max_abs) * 1.05  # 5% padding

    p <- ggplot2::ggplot(
        plot_df, 
        ggplot2::aes(x = Path, y = Estimate, ymin = Lower, ymax = Upper, color = Significant)
    ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::geom_pointrange(linewidth = 0.8, size = 0.5) +
        ggplot2::coord_flip(ylim = axis_lim) +
        ggplot2::scale_color_manual(values = color_map, guide = "none") +
        ggplot2::labs(
            title = paste("Path Coefficients (", type, ")", sep = ""),
            subtitle = if (type == "marginal") paste0("Average Marginal Effects with ", prob * 100, "% Bayesian Credibility Intervals") else paste0("Raw model parameters with ", prob * 100, "% Bayesian Credibility Intervals"),
            x = "",
            y = "Effect Size"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(family = "mono", size = 9),
            panel.grid.minor = ggplot2::element_blank()
        )

    return(p)
}
