#' @rdname pp_check
#'
#' @description
#' A wrapper around \code{bayesplot} posterior predictive check (PPC) functions
#' for \code{because} model objects.
#'
#' @param object A \code{because} fit object.
#' @param resp Character string; the response variable to check.
#'   If \code{NULL}, takes the first response.
#' @param type Character string; the type of PPC plot to generate.
#'   Supported: \code{"dens_overlay"}, \code{"hist"}, \code{"stat"}, \code{"stat_2d"},
#'   \code{"bars"}, \code{"bars_grouped"}, \code{"intervals"}, \code{"ecdf_overlay"},
#'   \code{"scatter_avg"}. Defaults to \code{"dens_overlay"}.
#' @param ndraws Integer; number of posterior draws to use. Defaults to 50.
#' @param trim Logical; if TRUE, zooms plot window to observed data range.
#'   Defaults to \code{FALSE} to allow inspection of predictive tail behavior.
#' @param re_formula Formula or \code{NA}; determines which random effects to include
#'   in the posterior predictions. See \code{\link{posterior_predict}} for details.
#' @param ... Additional arguments passed to \code{bayesplot} functions.
#'
#' @return A \code{ggplot} object produced by \code{bayesplot}.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(Y = rnorm(100), X = rnorm(100))
#' fit <- because(list(Y ~ X), data = df)
#'
#' # Default density overlay check
#' pp_check(fit)
#'
#' # Histogram check for a specific response
#' pp_check(fit, resp = "Y", type = "hist", ndraws = 30)
#'
#' # Test statistic check (e.g., mean)
#' pp_check(fit, type = "stat", stat = "mean")
#'
#' # Discrete bar plot check for count/binary models
#' pp_check(fit, type = "bars")
#' }
#'
#' @importFrom bayesplot ppc_dens_overlay ppc_hist ppc_stat ppc_stat_2d ppc_bars ppc_bars_grouped ppc_intervals ppc_ecdf_overlay ppc_scatter_avg
#' @importFrom ggplot2 ggplot
#' @export
pp_check.because <- function(object, resp = NULL, type = "dens_overlay", ndraws = 50, trim = FALSE, re_formula = NULL, ...) {
  # 1. Identify Response
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  # 2. Get Observed Data
  y_raw <- if (!is.null(object$original_data)) object$original_data[[resp]] else NULL
  if (is.null(y_raw)) {
     y_raw <- object$data[[resp]]
  }
  
  if (is.null(y_raw)) {
    stop(paste("Could not find observed data for response:", resp))
  }

  if (is.matrix(y_raw) && ncol(y_raw) == 1) {
    y_raw <- as.vector(y_raw)
  }
  
  # 3. Generate Posterior Predictions
  yrep <- posterior_predict(object, resp = resp, ndraws = ndraws, re_formula = re_formula)
  
  # Synchronous alignment and NA handling
  if (length(y_raw) == ncol(yrep)) {
    valid_idx <- which(!is.na(y_raw))
    if (length(valid_idx) == 0) {
      stop(paste("All observed values for response", resp, "are NA."))
    }
    y <- y_raw[valid_idx]
    yrep <- yrep[, valid_idx, drop = FALSE]
  } else if (length(y_raw) > ncol(yrep)) {
    # If yrep belongs to a higher level (e.g. 50 species) but y_raw is observation level (4500)
    h_info <- object$hierarchical_info
    target_lvl <- NULL
    if (!is.null(h_info)) {
      for (lvl in names(h_info$levels)) {
        if (resp %in% h_info$levels[[lvl]]) target_lvl <- lvl
      }
    }
    
    aligned <- FALSE
    # Check if original_data contains the level-specific dataframe
    if (!is.null(target_lvl) && is.list(object$original_data) && !is.data.frame(object$original_data)) {
      lvl_df <- object$original_data[[target_lvl]]
      if (is.data.frame(lvl_df) && resp %in% names(lvl_df)) {
        y_lvl <- lvl_df[[resp]]
        if (length(y_lvl) == ncol(yrep)) {
          valid_idx <- which(!is.na(y_lvl))
          y <- y_lvl[valid_idx]
          yrep <- yrep[, valid_idx, drop = FALSE]
          aligned <- TRUE
        }
      }
    }
    
    # Try finding unique entity values via linking variable
    if (!aligned && !is.null(target_lvl) && !is.null(h_info$link_vars[[target_lvl]])) {
      idx_var <- h_info$link_vars[[target_lvl]]
      links <- NULL
      if (idx_var %in% names(object$original_data)) {
        links <- object$original_data[[idx_var]]
      } else if (paste0(idx_var, "_idx") %in% names(object$data)) {
        links <- object$data[[paste0(idx_var, "_idx")]]
      }
      
      if (!is.null(links) && length(links) == length(y_raw)) {
        first_obs_idx <- !duplicated(links)
        y_subset <- y_raw[first_obs_idx]
        if (length(y_subset) == ncol(yrep)) {
          valid_idx <- which(!is.na(y_subset))
          y <- y_subset[valid_idx]
          yrep <- yrep[, valid_idx, drop = FALSE]
          aligned <- TRUE
        }
      }
    }
    
    if (!aligned) {
      stop(sprintf(
        "Cannot align observed data for '%s' (length %d) with posterior predictions (length %d). Observed data appears to be at a different hierarchical resolution.",
        resp, length(y_raw), ncol(yrep)
      ))
    }
  } else {
    stop(sprintf(
      "Posterior predictions has more columns (%d) than observed data (%d) for response '%s'.",
      ncol(yrep), length(y_raw), resp
    ))
  }

  # 4. Plot using bayesplot
  plot_out <- switch(type,
    "dens_overlay" = bayesplot::ppc_dens_overlay(y, yrep, ...),
    "hist"         = bayesplot::ppc_hist(y, yrep, ...),
    "stat"         = bayesplot::ppc_stat(y, yrep, ...),
    "stat_2d"      = bayesplot::ppc_stat_2d(y, yrep, ...),
    "bars"         = bayesplot::ppc_bars(y, yrep, ...),
    "bars_grouped" = bayesplot::ppc_bars_grouped(y, yrep, ...),
    "intervals"    = bayesplot::ppc_intervals(y, yrep, ...),
    "ecdf_overlay" = bayesplot::ppc_ecdf_overlay(y, yrep, ...),
    "scatter_avg"  = bayesplot::ppc_scatter_avg(y, yrep, ...),
    stop(paste("Unsupported PPC type:", type, ". Supported types include: 'dens_overlay', 'hist', 'stat', 'stat_2d', 'bars', 'bars_grouped', 'intervals', 'ecdf_overlay', 'scatter_avg'."))
  )
  
  # 5. Prettify with informative labels
  y_lab <- switch(type,
    "dens_overlay" = "Density",
    "hist"         = "Count",
    "stat"         = "Frequency",
    "stat_2d"      = "Stat 2",
    "bars"         = "Count",
    "bars_grouped" = "Count",
    "intervals"    = "Interval",
    "ecdf_overlay" = "ECDF",
    "scatter_avg"  = "Average yrep",
    "Value"
  )
  
  plot_out <- plot_out + 
    ggplot2::labs(
      title = paste("Posterior Predictive Check:", resp),
      subtitle = paste("Type:", type, "| draws:", nrow(yrep)),
      x = resp,
      y = y_lab
    ) +
    ggplot2::theme_minimal()

  # 6. Smart Zoom: Focus on observed data range (opt-in)
  if (trim) {
    y_range <- range(y, na.rm = TRUE)
    pad <- 0.05 * diff(y_range)
    if (pad == 0) pad <- 0.1
    
    plot_out <- plot_out + 
      ggplot2::coord_cartesian(xlim = c(y_range[1] - pad, y_range[2] + pad))
  }
  
  return(plot_out)
}

#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}
