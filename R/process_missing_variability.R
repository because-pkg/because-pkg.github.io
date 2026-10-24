`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Process Link Variables, Missing Data, Variability, and Latent Variables
#'
#' @noRd
process_missing_and_variability <- function(
  data,
  hierarchical_info,
  N,
  structures,
  family,
  family_obj,
  equations,
  fix_residual_variance,
  variability,
  latent,
  all_poly_terms,
  quiet,
  dsep,
  latent_method
) {
  # --- Map Link Variables (IDs) to Counts ---
  # Ensures JAGS finds loop bounds for grouping variables not listed in 'levels'
  if (!is.null(hierarchical_info$link_vars) && !is.null(hierarchical_info$data)) {
    for (lk_var in hierarchical_info$link_vars) {
      # Find which level (dataframe) contains this link variable
      # If multiple, find the one with the fewest rows (the level where it's defined)
      potential_lvls <- character(0)
      for (l_nm in names(hierarchical_info$data)) {
        if (lk_var %in% colnames(hierarchical_info$data[[l_nm]])) {
          potential_lvls <- c(potential_lvls, l_nm)
        }
      }

      if (length(potential_lvls) > 0) {
        # Pick the level with minimum rows
        r_counts <- vapply(potential_lvls, function(l) nrow(hierarchical_info$data[[l]]), numeric(1))
        best_lvl <- potential_lvls[which.min(r_counts)]
        best_n <- r_counts[best_lvl]

        # Generate names (Raw, Title, Upper)
        v_title <- paste0(toupper(substring(lk_var, 1, 1)), substring(lk_var, 2))
        pot_names <- unique(c(lk_var, v_title, toupper(lk_var)))

        for (p_nm in pot_names) {
          nn_name <- paste0("N_", p_nm)
          if (is.null(data[[nn_name]])) {
            data[[nn_name]] <- best_n
          }
        }
      }
    }
  }

  # Only add 'zeros' vector if using ZIP or ZINB (Poisson trick)
  # OR if we have structures (which often use dmnorm(zeros, ...))
  if (!is.null(N)) {
    has_zeros <- "zeros" %in% names(data)
    if (!has_zeros) {
      needs_zeros <- length(structures) > 0 ||
        any(vapply(
          names(family),
          function(v) {
            fam_name <- family[[v]]
            fam_obj_v <- get_family(fam_name)
            needs_zero_inflation_hook(fam_obj_v, v)
          },
          logical(1)
        ))

      if (needs_zeros) {
        data[["zeros"]] <- rep(0, N)
      }
    }
  }

  # ID2 is used for pairwise induced correlations (Wishart priors)
  data$ID2 <- diag(2)

  # Handle multinomial and ordinal data
  if (!is.null(family)) {
    # If family is provided but unnamed, try to auto-assign if there is only one response
    if (is.null(names(family))) {
      response_vars <- unique(vapply(
        equations,
        function(eq) as.character(all.vars(eq[[2]])[1]),
        character(1)
      ))
      if (length(family) == 1 && length(response_vars) == 1) {
        names(family) <- response_vars
        message(sprintf(
          "Auto-assigned family '%s' to response variable '%s'",
          family,
          response_vars
        ))
      } else if (length(family) == length(response_vars)) {
        warning(
          "Argument 'family' is unnamed. Please provide a named vector like c(Response = 'binomial'). Assuming defaults (Gaussian) for safety."
        )
      } else {
        warning(
          "Argument 'family' is unnamed and length does not match response variables. Ignoring."
        )
      }
    }

    # Auto-fix residual variance for non-Gaussian distributions if not specified
    for (var in names(family)) {
      dist <- family[[var]]
      if (dist %in% c("binomial", "multinomial", "ordinal")) {
        should_fix <- FALSE
        if (is.null(fix_residual_variance)) {
          should_fix <- TRUE
          fix_residual_variance <- c()
        } else if (
          is.numeric(fix_residual_variance) &&
            length(fix_residual_variance) == 1 &&
            is.null(names(fix_residual_variance))
        ) {
          should_fix <- FALSE
        } else if (!var %in% names(fix_residual_variance)) {
          should_fix <- TRUE
        }

        if (should_fix) {
          new_fix <- setNames(1, var)
          fix_residual_variance <- c(fix_residual_variance, new_fix)

          if (!quiet) {
            message(sprintf(
              "Note: Fixing residual variance of '%s' (%s) to 1 for identifiability.",
              var,
              dist
            ))
          }
        }
      }
    }

    # Identify categorical variables and their K levels
    cat_vars_metadata <- attr(data, "categorical_vars")

    for (var in names(family)) {
      if (family[[var]] %in% c("multinomial", "ordinal")) {
        if (!var %in% names(data)) {
          stop(paste(
            family[[var]],
            "variable",
            var,
            "not found in data."
          ))
        }

        # Determine K (number of levels)
        if (!is.null(cat_vars_metadata) && var %in% names(cat_vars_metadata)) {
          K <- length(cat_vars_metadata[[var]]$levels)
          val <- data[[var]]
          if (is.factor(val)) {
            data[[var]] <- as.integer(val)
          } else if (is.character(val)) {
            data[[var]] <- as.integer(factor(val, levels = cat_vars_metadata[[var]]$levels))
          }
        } else {
          val <- data[[var]]
          if (is.factor(val)) {
            K <- nlevels(val)
            data[[var]] <- as.integer(val)
          } else {
            val <- as.factor(val)
            K <- nlevels(val)
            data[[var]] <- as.integer(val)
          }
        }

        if (family[[var]] == "multinomial" && K < 3) {
          warning(paste(
            "Multinomial variable",
            var,
            "has fewer than 3 levels. Consider using binomial."
          ))
        }

        if (family[[var]] == "ordinal" && K < 3) {
          warning(paste(
            "Ordinal variable",
            var,
            "has fewer than 3 levels. Consider using binomial."
          ))
        }

        K_name <- paste0("K_", var)
        if (!K_name %in% names(data)) {
          data[[K_name]] <- K
          if (!quiet) {
            message(sprintf(
              "Auto-detected K_%s = %d from %s variable '%s'",
              var,
              K,
              family[[var]],
              var
            ))
          }
        }
      }
    }
  }

  if (!is.null(family)) {
    for (var_name in names(family)) {
      dist_type <- family[[var_name]]
      if (dist_type %in% c("zip", "zinb")) {
        zeros_name <- paste0("zeros_", var_name)
        if (is.null(data[[zeros_name]])) {
          data[[zeros_name]] <- rep(0, N)
        }
      }
    }
  }

  # Check for missing data
  all_vars <- unique(unlist(lapply(equations, function(eq) {
    c(all.vars(eq[[3]]), all.vars(eq[[2]]))
  })))

  cat_metadata <- attr(data, "categorical_vars")
  if (!is.null(cat_metadata)) {
    for (parent in names(cat_metadata)) {
      dummies <- cat_metadata[[parent]]$dummies
      if (any(dummies %in% all_vars)) {
        all_vars <- unique(c(all_vars, parent))
      }
    }
  }

  response_vars <- unique(vapply(
    equations,
    function(eq) as.character(all.vars(eq[[2]])[1]),
    character(1)
  ))
  predictor_only_vars <- setdiff(all_vars, response_vars)

  # Detect variables with missing data
  response_vars_with_na <- character(0)
  predictor_vars_with_na <- character(0)

  for (var in all_vars) {
    var_data <- if (var %in% names(data)) {
      data[[var]]
    } else if (is.list(data) && !is.data.frame(data)) {
      found_val <- NULL
      for (lvl_name in names(data)) {
        if (is.data.frame(data[[lvl_name]]) && var %in% colnames(data[[lvl_name]])) {
          found_val <- data[[lvl_name]][[var]]
          break
        }
      }
      found_val
    } else {
      NULL
    }

    if (!is.null(var_data)) {
      if (
        !is.matrix(var_data) && any(is.na(var_data)) && !all(is.na(var_data))
      ) {
        if (var %in% response_vars) {
          response_vars_with_na <- c(response_vars_with_na, var)
        } else {
          predictor_vars_with_na <- c(predictor_vars_with_na, var)
        }
      }
    }
  }

  # Handle predictor-only variables with missing data
  if (length(predictor_vars_with_na) > 0) {
    if (!quiet) {
      message(
        "Note: Detected missing data in predictor-only variables: ",
        paste(predictor_vars_with_na, collapse = ", "),
        "\nAutomatically adding intercept-only equations (e.g., ",
        predictor_vars_with_na[1],
        " ~ 1) to enable imputation."
      )
    }

    # Add intercept-only equations
    for (var in predictor_vars_with_na) {
      new_eq <- stats::as.formula(paste(var, "~ 1"))
      equations <- c(equations, list(new_eq))
    }

    # Treat them as responses now
    response_vars_with_na <- c(response_vars_with_na, predictor_vars_with_na)
  }

  # Auto-detect variability from data column names
  auto_variability <- list()

  for (var in all_vars) {
    if (!is.null(variability) && var %in% c(names(variability), variability)) {
      next
    }

    if (!is.null(family_obj)) {
      v_type <- get_variability_type_hook(family_obj, var)
      if (!is.null(v_type)) {
        auto_variability[[var]] <- v_type
        if (!quiet) {
          message(sprintf(
            "Extension-detected: '%s' is a specialized family -> using '%s' mode.",
            var,
            v_type
          ))
        }
        next
      }
    }

    se_name <- paste0(var, "_se")
    sd_name <- paste0(var, "_sd")

    if (se_name %in% names(data)) {
      auto_variability[[var]] <- "se"
      if (!quiet) {
        message(sprintf(
          "Auto-detected: '%s' has standard errors in '%s'",
          var,
          se_name
        ))
      }

      obs_name <- paste0(var, "_obs")
      if (var %in% names(data) && is.matrix(data[[var]])) {
        auto_variability[[var]] <- "reps"
        if (!quiet) {
          message(sprintf(
            "Auto-detected: '%s' has repeated measures (matrix format)",
            var
          ))
        }
      } else if (obs_name %in% names(data)) {
        auto_variability[[var]] <- "reps"
        if (!quiet) {
          message(sprintf(
            "Auto-detected: '%s' has repeated measures in '%s'",
            var,
            obs_name
          ))
        }
      }
    }
  }

  # Merge auto-detected with manual specification
  if (length(auto_variability) > 0) {
    if (is.null(variability)) {
      variability <- auto_variability
    } else {
      if (is.null(names(variability))) {
        variability <- setNames(rep("se", length(variability)), variability)
      }

      for (var in names(auto_variability)) {
        if (!var %in% names(variability)) {
          variability[[var]] <- auto_variability[[var]]
        }
      }
    }
  }

  # Handle variability data
  variability_list <- list()
  if (!is.null(variability)) {
    for (var_name in names(variability)) {
      var_spec <- variability[[var_name]]

      if (is.list(var_spec)) {
        type <- var_spec$type
        custom_se_col <- var_spec$se_col
        custom_obs_col <- var_spec$obs_col
        custom_mean_col <- var_spec$mean_col
      } else {
        type <- as.character(var_spec)
        custom_se_col <- NULL
        custom_obs_col <- NULL
        custom_mean_col <- NULL
      }

      if (!type %in% c("se", "reps")) {
        stop(paste(
          "Invalid variability type for",
          var_name,
          "- must be 'se' or 'reps', got:",
          type
        ))
      }

      variability_list[[var_name]] <- type

      if (type == "se") {
        se_col <- custom_se_col %||% paste0(var_name, "_se")
        mean_col <- custom_mean_col %||% paste0(var_name, "_mean")

        if (!se_col %in% names(data)) {
          stop(paste(
            "Variable",
            var_name,
            "specified as 'se' type but column",
            se_col,
            "not found in data."
          ))
        }

        if (!mean_col %in% names(data)) {
          if (var_name %in% names(data)) {
            data[[mean_col]] <- data[[var_name]]
            data[[var_name]] <- NULL
          } else {
            stop(paste(
              "Variable",
              var_name,
              "specified as 'se' type but neither",
              var_name,
              "nor",
              mean_col,
              "found in data."
            ))
          }
        } else {
          if (var_name %in% names(data)) data[[var_name]] <- NULL
        }

        if (se_col != paste0(var_name, "_se")) {
          data[[paste0(var_name, "_se")]] <- data[[se_col]]
        }
      } else if (type == "reps") {
        obs_col <- custom_obs_col %||% paste0(var_name, "_obs")
        nrep_name <- paste0("N_reps_", var_name)

        if (!obs_col %in% names(data)) {
          if (
            var_name %in%
              names(data) &&
              (is.matrix(data[[var_name]]) || is.data.frame(data[[var_name]]))
          ) {
            data[[obs_col]] <- as.matrix(data[[var_name]])
            data[[var_name]] <- NULL
          } else if (
            var_name %in%
              names(data) &&
              is.list(data[[var_name]]) &&
              !is.data.frame(data[[var_name]])
          ) {
            elem_valid <- all(sapply(data[[var_name]], function(x) {
              is.matrix(x) || is.data.frame(x)
            }))
            if (!elem_valid) {
              stop(
                "Elements of '",
                var_name,
                "' list must be matrices or data frames."
              )
            }

            tryCatch(
              {
                mat_list <- lapply(data[[var_name]], as.matrix)
                arr_3d <- simplify2array(mat_list)
                data[[obs_col]] <- arr_3d
                data[[var_name]] <- NULL

                if (!quiet) {
                  message(
                    "Converted list of matrices '",
                    var_name,
                    "' to 3D array (",
                    paste(dim(arr_3d), collapse = "x"),
                    ")."
                  )
                }
              },
              error = function(e) {
                stop(
                  "Failed to convert list '",
                  var_name,
                  "' to 3D array. Ensure all matrices have identical dimensions. Error: ",
                  e$message
                )
              }
            )
          } else {
            stop(paste(
              "Variable",
              var_name,
              "specified as 'reps' type but column",
              obs_col,
              "(as matrix) not found in data."
            ))
          }
        }

        if (obs_col != paste0(var_name, "_obs")) {
          data[[paste0(var_name, "_obs")]] <- data[[obs_col]]
        }

        if (!nrep_name %in% names(data)) {
          mat <- data[[paste0(var_name, "_obs")]]
          n_reps <- apply(mat, 1, function(x) sum(!is.na(x)))
          data[[nrep_name]] <- n_reps

          compact_mat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
          for (i in seq_len(nrow(mat))) {
            vals <- mat[i, !is.na(mat[i, ])]
            if (length(vals) > 0) {
              compact_mat[i, seq_along(vals)] <- vals
            }
          }
          data[[paste0(var_name, "_obs")]] <- compact_mat
        }
      }

      if (var_name %in% names(data)) data[[var_name]] <- NULL
    }
  }

  # Auto-detect latent variables: variables in equations but not in data
  if (is.null(latent)) {
    vars_in_equations <- unique(unlist(lapply(equations, all.vars)))
    vars_in_data <- names(data)

    potential_latents <- setdiff(vars_in_equations, vars_in_data)

    potential_latents <- dsep_potential_latent_hook(
      family_obj,
      potential_latents
    )

    if (!is.null(all_poly_terms)) {
      poly_internal_names <- sapply(all_poly_terms, function(x) {
        x$internal_name
      })
      potential_latents <- setdiff(potential_latents, poly_internal_names)
    }

    cat_metadata <- attr(data, "categorical_vars")
    if (!is.null(cat_metadata)) {
      potential_latents <- setdiff(potential_latents, names(cat_metadata))
    }

    if (length(potential_latents) > 0) {
      latent <- potential_latents

      if (!quiet) {
        msg <- paste0(
          "Auto-detected latent variable(s): ",
          paste(latent, collapse = ", "),
          "\n(Variables in equations but not in data will be treated as latent.)"
        )
        if (dsep) {
          msg <- paste0(msg, "\nGenerating m-separation tests for MAG...")
        }
        message(msg)
      }
    }
  }

  if (!quiet && dsep) {
    if (!is.null(latent)) {
      message("Generating m-separation tests (MAG with latent variables)...")
    } else {
      message("Generating d-separation tests...")
    }
  }

  list(
    data                  = data,
    family                = family,
    fix_residual_variance = fix_residual_variance,
    equations             = equations,
    variability           = variability,
    variability_list      = variability_list,
    latent                = latent,
    response_vars_with_na = response_vars_with_na
  )
}
