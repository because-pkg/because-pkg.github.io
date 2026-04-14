#' Posterior Predictive Checks for Because Fits
#'
#' @description
#' A wrapper around \code{bayesplot::ppc_dens_overlay} and other PPC functions
#' for \code{because} model objects.
#'
#' @param object A \code{because} fit object.
#' @param resp Character string; the response variable to check.
#'   If \code{NULL}, takes the first response.
#' @param type Character string; the type of PPC plot to generate.
#'   Supported: \code{"dens_overlay"}, \code{"hist"}, \code{"stat"}.
#' @param ndraws Integer; number of posterior draws to use. Defaults to 50.
#' @param ... Additional arguments passed to \code{bayesplot} functions.
#'
#' @return A \code{ggplot} object produced by \code{bayesplot}.
#'
#' @importFrom bayesplot ppc_dens_overlay ppc_hist ppc_stat
#' @importFrom ggplot2 ggplot
#' @export
pp_check.because <- function(object, resp = NULL, type = "dens_overlay", ndraws = 50, ...) {
  # 1. Identify Response
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  # 2. Get Observed Data
  y <- object$original_data[[resp]]
  if (is.null(y)) {
     y <- object$data[[resp]]
  }
  
  if (is.null(y)) {
    stop(paste("Could not find observed data for response:", resp))
  }
  
  # Remove NAs for plotting
  y <- na.omit(y)
  
  # 3. Generate Posterior Predictions
  yrep <- posterior_predict(object, resp = resp, ndraws = ndraws)
  
  # Ensure dimensions match (in case data was filtered during fit)
  if (ncol(yrep) != length(y)) {
    # If using processed data, it might have been filtered for NAs
    # We should match based on names if available
    message("Note: y and yrep dimensions mismatch. Attempting to align...")
    if (length(y) > ncol(yrep)) {
        y <- y[seq_len(ncol(yrep))]
    } else {
        yrep <- yrep[, seq_along(y), drop=FALSE]
    }
  }

  # 4. Plot using bayesplot
  plot_out <- switch(type,
    "dens_overlay" = bayesplot::ppc_dens_overlay(y, yrep, ...),
    "hist"         = bayesplot::ppc_hist(y, yrep, ...),
    "stat"         = bayesplot::ppc_stat(y, yrep, ...),
    stop(paste("Unknown PPC type:", type))
  )
  
  # 5. Prettify with informative labels
  y_lab <- switch(type,
    "dens_overlay" = "Density",
    "hist"         = "Count",
    "stat"         = "Frequency",
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
  
  return(plot_out)
}

#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}
