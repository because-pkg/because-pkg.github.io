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
pp_check.because <- function(object, resp = NULL, type = "dens_overlay", ndraws = 50, trim = TRUE, ...) {
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
  
  # Ensure dimensions match
  if (ncol(yrep) != length(y)) {
    # If yrep belongs to a higher level (e.g. 50 species) but y is observation level (4500)
    # we need to subset y to match.
    if (length(y) > ncol(yrep)) {
        h_info <- object$hierarchical_info
        target_lvl <- NULL
        if (!is.null(h_info)) {
            for (lvl in names(h_info$levels)) if (resp %in% h_info$levels[[lvl]]) target_lvl <- lvl
        }
        
        # If we found a level, try to find the index to subset y
        subset_done <- FALSE
        if (!is.null(target_lvl)) {
            idx_var <- h_info$link_vars[[target_lvl]]
            
            # Find the linking column (search original_data list or flat data)
            links <- NULL
            if (!is.null(idx_var)) {
                if (idx_var %in% names(object$original_data)) {
                    links <- object$original_data[[idx_var]]
                } else if (is.list(object$original_data)) {
                    # Search inside sub-dataframes (usually 'obs' has the links)
                    for (df_name in names(object$original_data)) {
                        if (is.data.frame(object$original_data[[df_name]]) && idx_var %in% names(object$original_data[[df_name]])) {
                            links <- object$original_data[[df_name]][[idx_var]]
                            break
                        }
                    }
                }
                
                if (is.null(links) && paste0(idx_var, "_idx") %in% names(object$data)) {
                    links <- object$data[[paste0(idx_var, "_idx")]]
                }
            }
            
            if (!is.null(links)) {
                # Ensure links and y have same length
                if (length(links) == length(y)) {
                    # Take first observation per entity
                    first_obs_idx <- !duplicated(links)
                    y_subset <- y[first_obs_idx]
                    if (length(y_subset) == ncol(yrep)) {
                        y <- y_subset
                        subset_done <- TRUE
                    }
                }
            }
        }
        
        if (!subset_done) {
            # Fallback to simple truncation
            message("Note: y and yrep dimensions mismatch. Attempting to align via truncation...")
            y <- y[seq_len(ncol(yrep))]
        }
    } else {
        # yrep is larger? Should not happen with within-sample PPC
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

  # 6. Smart Zoom: Focus on observed data range
  if (trim) {
    y_range <- range(y, na.rm = TRUE)
    # Expand slightly
    pad <- 0.05 * diff(y_range)
    if (pad == 0) pad <- 0.1 # Constant if no variance
    
    plot_out <- plot_out + 
      ggplot2::coord_cartesian(xlim = c(y_range[1] - pad, y_range[2] + pad))
  }
  
  return(plot_out)
}

#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}
