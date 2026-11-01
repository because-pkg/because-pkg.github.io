#' Check for unstandardized continuous variables with large variance
#'
#' Inspects continuous (Gaussian) response variables in the model equations.
#' If any have large empirical standard deviations (SD > 20) and no custom
#' prior was provided for that variable, issues an informative warning
#' alerting the user that default priors assume standardized variables.
#'
#' @param data Model data (data.frame or list)
#' @param equations List of model equations (formulas)
#' @param family Named vector or list of distribution families
#' @param priors List of custom priors provided by the user (or NULL)
#' @param quiet Logical; if TRUE, suppress warnings
#' @keywords internal
check_prior_scale <- function(data, equations, family = NULL, priors = NULL, quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))
  
  # Extract response variables that are modeled with Gaussian distribution
  gaussian_responses <- character(0)
  for (eq in equations) {
    f_eq <- tryCatch(stats::formula(eq), error = function(e) NULL)
    if (is.null(f_eq) || length(f_eq) < 2) next
    resp <- as.character(f_eq[[2]])[1]
    fam <- if (!is.null(family[[resp]])) tolower(family[[resp]]) else "gaussian"
    if (fam %in% c("gaussian", "normal")) {
      gaussian_responses <- c(gaussian_responses, resp)
    }
  }
  
  if (length(gaussian_responses) == 0) return(invisible(NULL))
  
  # Helper to retrieve vector from data (handles data.frame or list of data.frames)
  get_vec <- function(vname) {
    if (is.data.frame(data) && vname %in% names(data)) {
      return(data[[vname]])
    } else if (is.list(data)) {
      if (vname %in% names(data) && !is.data.frame(data[[vname]])) {
        return(data[[vname]])
      }
      for (item in data) {
        if (is.data.frame(item) && vname %in% names(item)) {
          return(item[[vname]])
        }
      }
    }
    return(NULL)
  }
  
  unstandardized <- list()
  for (resp in unique(gaussian_responses)) {
    # Check if a custom prior was already supplied for this response's variance/precision
    if (!is.null(priors)) {
      prior_names <- names(priors)
      if (any(grepl(paste0("(sigma|tau).*", resp), prior_names))) {
        next
      }
    }
    
    vec <- get_vec(resp)
    if (is.numeric(vec)) {
      s_val <- stats::sd(vec, na.rm = TRUE)
      if (!is.na(s_val) && s_val > 20) {
        unstandardized[[resp]] <- s_val
      }
    }
  }
  
  if (length(unstandardized) > 0) {
    msg_vars <- paste(vapply(names(unstandardized), function(v) {
      sprintf("'%s' (SD = %.1f)", v, unstandardized[[v]])
    }, character(1)), collapse = ", ")
    
    warning(
      sprintf(
        "Continuous variable(s) with large variance detected: %s. Default priors in because (e.g., sigma ~ Uniform(0, 100), alpha/beta ~ Normal(0, 100)) assume standardized variables (mean = 0, SD = 1). Please consider standardizing continuous variables or supplying scale-appropriate custom priors via the 'priors' argument to prevent unintended prior truncation.",
        msg_vars
      ),
      call. = FALSE
    )
  }
  
  invisible(NULL)
}
