#' Plot Posterior Distributions for Model Comparison
#'
#' Visualizes posterior density distributions for selected parameters, allowing easy comparison
#' between multiple models (e.g., Default vs. Custom priors) or inspection of a single model.
#'
#' @param models A `because` object or a named list of `because` objects.
#'   If a list is provided, the names are used for the legend.
#'   Example: `list("Default" = fit1, "Custom" = fit2)`.
#' @param parameter Character string. A regular expression or exact name of the parameter(s) to plot.
#'   Example: `"beta"`, `"^alpha"`, `"beta_Y_X"`. Plots all matching parameters.
#' @param col Optional vector of colors. Defaults to a standard palette (black, blue, red, etc.).
#' @param lwd Line width (default = 2).
#' @param legend_pos Legend position (default = "topleft"). Set to `NULL` to suppress.
#' @param density_args List of additional arguments passed to `stats::density()`.
#'   Useful for handling boundary corrections (e.g., `list(from=0)` for truncated priors).
#' @param ... Additional arguments passed to `plot()`.
#'
#' @return Invisible NULL. Produces a plot.
#' @examples
#' \dontrun{
#' fit <- because(list(Y ~ X), data = my_data)
#' plot_posterior(fit, parameter = "beta_Y_X")
#' }
#'
#' @export
#' @importFrom graphics plot lines legend par grid axis layout
#' @importFrom stats density
#' @importFrom grDevices rainbow palette
plot_posterior <- function(
    models,
    parameter,
    col = NULL,
    lwd = 2,
    legend_pos = "topleft",
    density_args = list(),
    ...
) {
    # Standardize input to list
    if (inherits(models, "because")) {
        model_list <- list("Model" = models)
    } else if (is.list(models) && all(sapply(models, inherits, "because"))) {
        model_list <- models
        if (is.null(names(model_list))) {
            names(model_list) <- paste("Model", seq_along(model_list))
        }
    } else {
        stop(
            "Argument 'models' must be a 'because' object or a list of 'because' objects."
        )
    }

    # Colors
    if (is.null(col)) {
        # Simple palette
        def_cols <- c("black", "blue", "red", "green4", "purple", "orange")
        col <- rep(def_cols, length.out = length(model_list))
    } else {
        col <- rep(col, length.out = length(model_list))
    }

    # Identify parameters to plot (based on first model)
    # Intersection of parameters across all models would be safer, but first is reasonable anchor
    ref_model <- model_list[[1]]
    if (is.null(ref_model$samples)) {
        stop("Model(s) do not contain samples. Did you run with n.iter=0?")
    }

    all_params <- colnames(as.matrix(ref_model$samples[[1]]))
    # Match parameter regex
    targets <- grep(parameter, all_params, value = TRUE)

    if (length(targets) == 0) {
        stop(sprintf(
            "No parameters matched '%s'. Available params:\n%s",
            parameter,
            paste(head(all_params, 10), collapse = ", ")
        ))
    }

    # Setup layout if multiple parameters
    n_targets <- length(targets)
    if (n_targets > 1) {
        # Simple grid logic
        ncols <- ceiling(sqrt(n_targets))
        nrows <- ceiling(n_targets / ncols)
        old_par <- par(mfrow = c(nrows, ncols), mar = c(4, 4, 2, 1))
        on.exit(par(old_par))
    }

    # Loop through parameters
    for (p in targets) {
        # Collect density objects and ranges
        densities <- list()
        x_ranges <- c()
        y_ranges <- c()

        for (i in seq_along(model_list)) {
            m <- model_list[[i]]

            # Check if param exists in this model
            m_params <- colnames(as.matrix(m$samples[[1]]))
            if (!p %in% m_params) {
                warning(sprintf(
                    "Parameter '%s' not found in model '%s'. Skipping.",
                    p,
                    names(model_list)[i]
                ))
                densities[[i]] <- NULL
                next
            }

            # Extract samples
            samps <- as.matrix(m$samples)[, p]

            # Compute density with optional arguments
            # Check for per-model arguments
            current_args <- list()

            # Check if density_args implies named list of lists matching model names
            is_named_list <- !is.null(names(density_args)) &&
                any(names(density_args) %in% names(model_list))

            if (is_named_list) {
                # Look for args for this specific model name
                model_name <- names(model_list)[i]
                if (model_name %in% names(density_args)) {
                    current_args <- density_args[[model_name]]
                }
            } else {
                # Apply globally
                current_args <- density_args
            }

            args <- c(list(x = samps), current_args)
            d <- do.call(density, args)

            # Store computed density AND args used (for visualization decisions)
            densities[[i]] <- list(d = d, args = current_args)
            x_ranges <- c(x_ranges, d$x)
            y_ranges <- c(y_ranges, d$y)
        }

        # Plotting canvas
        if (length(x_ranges) == 0) {
            next
        } # Skip if no data found

        x_lim <- range(x_ranges)
        y_lim <- range(y_ranges)

        # Base Plot
        plot(
            NULL,
            xlim = x_lim,
            ylim = y_lim,
            main = p,
            xlab = "Estimate",
            ylab = "Density",
            ...
        )
        grid()

        # Add lines
        valid_indices <- which(!sapply(densities, is.null))
        for (i in valid_indices) {
            d_obj <- densities[[i]]$d
            args <- densities[[i]]$args

            x_vals <- d_obj$x
            y_vals <- d_obj$y

            # Add vertical drops if truncation arguments were used
            if (!is.null(args[["from"]])) {
                # Prepend (x[1], 0)
                x_vals <- c(x_vals[1], x_vals)
                y_vals <- c(0, y_vals)
            }
            if (!is.null(args[["to"]])) {
                # Append (x[n], 0)
                n <- length(x_vals)
                x_vals <- c(x_vals, x_vals[n])
                y_vals <- c(y_vals, 0)
            }

            lines(x_vals, y_vals, col = col[i], lwd = lwd)
        }

        # Legend (only on first plot if multiple)
        # Or maybe on all? Let's put on all for clarity, user can suppress with legend_pos=NULL
        if (!is.null(legend_pos) && length(model_list) > 1) {
            legend(
                legend_pos,
                legend = names(model_list),
                col = col,
                lwd = lwd,
                bty = "n",
                cex = 0.8
            )
        }
    }
}
