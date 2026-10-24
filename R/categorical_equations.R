#' Expand Categorical Variables in Structural Equations
#'
#' Automatically expands categorical variables in model formulas into dummy indicator
#' variables, supporting interaction terms and deterministic identity definitions
#' for imputed factors.
#'
#' @keywords internal
expand_categorical_equations <- function(
  equations,
  data,
  response_vars_with_na = NULL,
  quiet = FALSE
) {
  if (is.null(attr(data, "categorical_vars"))) {
    return(list(equations = equations, data = data))
  }

  categorical_vars <- attr(data, "categorical_vars")
  new_equations <- list()

  # Use for loop instead of lapply to safely modify data and collect equations
  for (idx in seq_along(equations)) {
    eq <- equations[[idx]]
    vars <- all.vars(eq)

    # Check if any predictors are categorical
    for (var in vars) {
      if (var %in% names(categorical_vars)) {
        # Check if var is the response (LHS)
        lhs_var <- all.vars(eq[[2]])
        if (var %in% lhs_var) {
          # Skip expansion if it's the response
          next
        }

        # Get dummy variable names
        levels <- categorical_vars[[var]]$levels
        dummies <- categorical_vars[[var]]$dummies

        # If the parent variable is being imputed (is in response_vars_with_na),
        # we MUST define the dummies deterministically in JAGS to link them.
        if (!is.null(response_vars_with_na) && var %in% response_vars_with_na) {
          for (k in 2:length(levels)) {
            dummy_name <- paste0(var, "_", levels[k])

            # Create deterministic equation: dummy ~ I(var == k)
            det_eq_str <- sprintf("%s ~ I(%s == %d)", dummy_name, var, k)
            det_eq <- stats::as.formula(det_eq_str)

            # Only add if not already present
            eq_exists <- any(sapply(c(equations, new_equations), function(e) {
              deparse(e) == deparse(det_eq)
            }))

            if (!eq_exists) {
              new_equations <- c(new_equations, list(det_eq))
              # Also remove the dummy from the data list so JAGS uses the definition
              if (dummy_name %in% names(data)) {
                data[[dummy_name]] <- NULL
              }
            }
          }
        }

        # Skip substitution if we are inside a deterministic identity definition
        # for an imputed factor (e.g., Rate_2 ~ I(Rate == 2)).
        # This prevents circular expansion of 'Rate' into its own dummies.
        is_identity_mapping <- length(vars) == 2 && grepl(sprintf("I\\(%s == [0-9]+\\)", var), deparse(eq))
        if (is_identity_mapping) {
          next
        }

        # Convert formula to character for manipulation
        eq_str <- paste(deparse(eq), collapse = " ")

        # Replace categorical variable with its dummies (wrapped in parentheses for proper interaction expansion like A*B)
        pattern <- paste0("\\b", var, "\\b")
        replacement <- paste0("(", paste(dummies, collapse = " + "), ")")
        eq_str <- gsub(pattern, replacement, eq_str)

        # Convert back to formula
        eq <- stats::as.formula(eq_str)
        equations[[idx]] <- eq

        if (!quiet) {
          message(sprintf(
            "Expanded '%s' to: %s",
            var,
            paste(dummies, collapse = ", ")
          ))
        }
      }
    }
  }

  # Add the new deterministic equations
  equations <- c(equations, new_equations)

  return(list(equations = equations, data = data))
}
