#' Process Latent Variable Structural Equations for MAG Logic
#'
#' Marginalizes latent variables from structural equations, computes induced
#' correlations if needed, and injects intercept-only models for disconnected variables.
#'
#' @keywords internal
process_latent_mag_equations <- function(
  latent,
  latent_method,
  equations,
  random_terms,
  hierarchical_info,
  family,
  quiet = FALSE,
  induced_cors = NULL,
  dsep = FALSE
) {
  if (is.null(latent)) {
    return(list(equations = equations, induced_cors = induced_cors))
  }

  latent_method <- match.arg(latent_method, c("correlations", "explicit"))

  # Force MAG approach when doing d-separation testing
  if (dsep && latent_method == "explicit") {
    if (!quiet) {
      message(
        "Note: d-separation testing with latent variables requires MAG approach. ",
        "Using latent_method = 'correlations'."
      )
    }
    latent_method <- "correlations"
  }

  if (latent_method == "correlations") {
    # MAG approach: marginalize latents, use induced correlations
    # If not already computed by dsep, compute now
    if (is.null(induced_cors)) {
      dsep_result <- because_dsep(
        equations,
        latent = latent,
        random_terms = random_terms,
        hierarchical_info = hierarchical_info,
        family = if (!is.null(family)) as.list(family) else NULL,
        quiet = !dsep
      )
      induced_cors <- dsep_result$correlations
    }

    # Filter out equations involving latent variables
    equations <- lapply(equations, function(eq) {
      vars <- all.vars(eq)
      if (any(vars %in% latent)) {
        # Remove latent variables from the formula
        rhs <- labels(stats::terms(eq))
        keep_terms <- rhs[!rhs %in% latent]

        if (length(keep_terms) == 0) {
          # Becomes intercept-only model
          new_eq <- stats::as.formula(paste(as.character(eq)[2], "~ 1"))
        } else {
          # Keep observed predictors
          new_eq <- stats::as.formula(paste(
            as.character(eq)[2],
            "~",
            paste(keep_terms, collapse = " + ")
          ))
        }
        return(new_eq)
      }
      return(eq)
    })

    if (!quiet) {
      message(
        "Using MAG approach: marginalized latent variables from structural equations."
      )
    }

    # Check if variables with induced correlations need intercept-only models
    # This happens when a variable is involved in an induced correlation
    # but is not a response variable in any remaining equation
    if (length(induced_cors) > 0) {
      vars_with_correlations <- unique(unlist(induced_cors))
      response_vars <- sapply(equations, function(eq) as.character(eq)[2])
      vars_needing_intercept <- setdiff(vars_with_correlations, response_vars)

      if (length(vars_needing_intercept) > 0) {
        intercept_equations <- lapply(vars_needing_intercept, function(v) {
          stats::as.formula(paste(v, "~ 1"))
        })
        equations <- c(equations, intercept_equations)

        if (!quiet) {
          message(
            "Created intercept-only models for ",
            length(vars_needing_intercept),
            " variable(s) with induced correlations: ",
            paste(vars_needing_intercept, collapse = ", ")
          )
        }
      }
    }

    if (!quiet && length(induced_cors) > 0) {
      message(
        "Estimating ",
        length(induced_cors),
        " induced correlation(s) from latent variable(s)"
      )
    }
  } else {
    # Explicit approach: keep all equations, don't use induced correlations
    induced_cors <- NULL

    if (!quiet) {
      message("Using explicit latent variable modeling")
    }
  }

  return(list(equations = equations, induced_cors = induced_cors))
}
