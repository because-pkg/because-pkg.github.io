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

#' @rdname plot_rope
#'
#' @description
#' @param object A `because` object.
#' @param rope A numeric vector of length 2 specifying the ROPE limits. Default is `c(-0.1, 0.1)`.
#' @param parameters Optional character vector to filter which parameters to plot.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `ggplot` object.
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_manual geom_vline labs theme_minimal theme element_text after_stat
plot_rope.because <- function(object, rope = c(-0.1, 0.1), parameters = NULL, ...) {
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
                        
                        vals <- samps[, param_name]
                        tmp <- data.frame(
                            Label = clean_label,
                            Value = vals
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
            target_params <- intersect(target_params, parameters)
        }
        
        if (length(target_params) == 0) {
            stop("No valid parameters found to plot.")
        }
        
        all_samps <- as.matrix(object$samples)
        
        for (p in target_params) {
            if (p %in% colnames(all_samps)) {
                vals <- all_samps[, p]
                tmp <- data.frame(
                    Label = p,
                    Value = vals
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
    
    rope_min <- rope[1]
    rope_max <- rope[2]
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, y = Label)) +
        ggridges::stat_density_ridges(
            geom = "density_ridges_gradient",
            ggplot2::aes(fill = ggplot2::after_stat(x >= rope_min & x <= rope_max)),
            color = "black",
            scale = 0.95,
            rel_min_height = 0.005,
            show.legend = TRUE
        ) +
        ggplot2::scale_fill_manual(
            name = "Region",
            values = c("TRUE" = "#89C5C5", "FALSE" = "#F8766D"),
            labels = c(
                "TRUE" = paste0("Inside ROPE [", rope_min, ", ", rope_max, "]"), 
                "FALSE" = "Outside ROPE"
            ),
            breaks = c("TRUE", "FALSE")
        ) +
        ggplot2::geom_vline(xintercept = rope, linetype = "dashed", color = "darkgray", linewidth = 1) +
        ggplot2::labs(
            title = if(is_dsep) "d-separation Tests: ROPE Overlap" else "Posterior Distributions: ROPE Overlap",
            subtitle = "Density curves stacked vertically",
            x = "Effect Size",
            y = ""
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = "bottom"
        )
        
    return(p)
}
