#' plot_coef
#'
#' Creates a caterpillar plot (point and whisker) of the path coefficients
#' from a `because` model. This provides a complementary view to the DAG,
#' showing the exact magnitudes and uncertainties of causal effects.
#'
#' @param object A `because` object.
#' @param type Character; either `"marginal"` (default) or `"raw"`.
#'   - `"marginal"`: Shows Average Marginal Effects (AME) on the response scale.
#'   - `"raw"`: Shows the raw structural parameters (betas/rhos).
#' @param multinomial_probabilities Logical; if `TRUE` (default), expands 
#'   multinomial predictors into a "bundle" of category-specific effects.
#' @param color_scheme Character; color scheme for significance. Options:
#'   - `"directional"` (default): Blue for positive significant, Red for negative, Grey otherwise.
#'   - `"sig_only"`: Black for significant, Grey for non-significant.
#'   - `"monochrome"`: All black.
#' @param ... Additional arguments.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_hline coord_flip labs theme_minimal scale_color_manual
plot_coef.because <- function(
    object,
    type = "marginal",
    multinomial_probabilities = TRUE,
    color_scheme = "directional",
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
        # Use a higher sample size for the coefficient plot than for the DAG (better precision)
        me_table <- marginal_effects(
            object, 
            samples = 250, 
            multinomial_probabilities = multinomial_probabilities
        )
    } else {
        if (is.null(object$summary)) {
            stats <- summary(object$samples)$statistics
            quantiles <- summary(object$samples)$quantiles
        } else {
            stats <- object$summary$statistics
            quantiles <- object$summary$quantiles
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
         target_params <- target_params[!grepl("(_det_|_times_)", target_params)]

         for (pname in target_params) {
             est  <- stats[pname, "Mean"]
             low  <- quantiles[pname, "2.5%"]
             upp  <- quantiles[pname, "97.5%"]

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
        "pos" = "#1f78b4",     # Blue
        "neg" = "#e31a1c",     # Red
        "sig" = "black",
        "non-sig" = "#999999"  # Grey
    )

    p <- ggplot2::ggplot(
        plot_df, 
        ggplot2::aes(x = Path, y = Estimate, ymin = Lower, ymax = Upper, color = Significant)
    ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        ggplot2::geom_pointrange(linewidth = 0.8, size = 0.5) +
        ggplot2::coord_flip() +
        ggplot2::scale_color_manual(values = color_map, guide = "none") +
        ggplot2::labs(
            title = paste("Path Coefficients (", type, ")", sep = ""),
            subtitle = if (type == "marginal") "Average Marginal Effects with 95% Bayesian Credibility Intervals" else "Raw model parameters",
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
