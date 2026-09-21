#' Prepare random effects data structures for because models
#'
#' Expands random effects terms into the data list, creating index vectors,
#' grouping matrices, and associated loop bounds needed by JAGS/NIMBLE.
#' @keywords internal
prepare_random_effects_data <- function(
  data, random_terms, equations, hierarchical_info, is_hierarchical,
  levels, family, quiet,
  variability = NULL, id_col = NULL, all_poly_terms = NULL,
  latent = NULL, structure = NULL
) {
  # --- Random Effects Data Prep (Post-Assembly) ---
  # Create structures for JAGS using the assembled data
  if (length(random_terms) > 0) {
    rand_structs <- create_group_structures(data, random_terms)

    if (!quiet) {}
    rand_structs <- create_group_structures(data, random_terms)
    random_structures <- rand_structs$structures
    random_data_updates <- rand_structs$data_updates
  }

  if (is.data.frame(data) || (is.list(data) && !is.data.frame(data))) {
    # Check for long format data requiring matrix conversion
    # If any variability specified as 'reps', we attempt to auto-format using because_format_data
    has_reps <- any(grepl("reps", variability))

    if (has_reps) {
      if (is.null(structure)) {
        # Cannot auto-format without structure to determine species order
        warning(
          "Variability 'reps' specified but no structure provided. Automatic formatting requires a structure to order species rows. Assuming data is already aggregated or user handles index mapping."
        )
        return(data)
      } else if (!is.null(id_col) && id_col %in% names(data)) {
        # Extension Hook: Extract appropriate tree from structure object
        use_tree <- get_tree_hook(structure)

        formatted_list <- because_format_data(
          data,
          species_col = id_col,
          tree = use_tree
        )

        # Now rename the variables that are 'reps' to include '_obs' suffix
        # And keep others as is
        reps_vars <- names(variability)[variability == "reps"]

        final_data_list <- list()
        for (nm in names(formatted_list)) {
          if (nm %in% reps_vars) {
            # This is a matrix of replicates, rename to _obs
            final_data_list[[paste0(nm, "_obs")]] <- formatted_list[[nm]]
          } else {
            # This is a regular variable (vector or matrix depending on because_format_data logic)
            final_data_list[[nm]] <- formatted_list[[nm]]
          }
        }

        # Update data to be this list
        data <- final_data_list

        if (!quiet) {
          message(
            "  Formatted ",
            length(names(formatted_list)),
            " variables as replicate matrices."
          )
        }

        # Check for missing variables that were expected to be formatted
        # Filter variability/family to variables actually in current equations
        # to avoid 'missing reps' errors for irrelevant variables in sub-models (e.g. d-sep tests)
        all_vars <- unique(c(
          names(equations),
          unlist(lapply(equations, all.vars))
        ))
        # Include versions with p_ removed for variability/dist matching
        clean_vars <- sub("^p_", "", all_vars)
        relevant_vars <- unique(c(all_vars, clean_vars))

        missing_reps <- setdiff(reps_vars, names(formatted_list))
        # Only error if the missing variable is actually relevant to our equations
        missing_reps <- intersect(missing_reps, relevant_vars)

        if (length(missing_reps) > 0) {
          stop(paste(
            "The following variables were identified for 'reps' processing (from equations) but were NOT found in the data:",
            paste(missing_reps, collapse = ", "),
            "\nPlease check your column names."
          ))
        }
      } else {
        warning(
          "Variability 'reps' specified but 'id_col' missing or not in data. Cannot auto-format long data."
        )
      }
    }
  }

  # row_ids was already captured above (immediately after hierarchical_info was
  # built from the raw data). For the flat path it is re-initialised below.
  # row_ids <- NULL  # <-- original position; now handled earlier.

  if (
    (is.data.frame(data) || (is.list(data) && !is.data.frame(data))) &&
      !is_hierarchical
  ) {
    # Extract all variable names from fixed equations
    eq_vars <- unique(unlist(lapply(equations, all.vars)))

    # Add variables from random terms (grouping factors)
    if (length(random_terms) > 0) {
      random_vars <- unique(vapply(
        random_terms,
        function(x) x$group,
        character(1)
      ))
      eq_vars <- unique(c(eq_vars, random_vars))
    }

    # Check which variables are in the data frame
    # Check which variables are in the data frame
    # We must include "N" if it exists, as it's needed for optimized loops
    cols_to_keep <- unique(c(eq_vars, "N"))
    if (!is.null(variability) && !is.character(variability)) {
      cols_to_keep <- c(cols_to_keep, names(variability))
    }

    available_vars <- intersect(
      cols_to_keep,
      (if (is.null(names(data))) character(0) else names(data))
    )
    missing_vars <- setdiff(
      eq_vars,
      (if (is.null(names(data))) character(0) else names(data))
    )

    # Some "missing" vars might be latent - that's OK
    if (!is.null(latent)) {
      missing_vars <- setdiff(missing_vars, latent)
    }

    if (length(missing_vars) > 0 && length(available_vars) == 0) {
      stop(
        "None of the variables in equations found in data frame. ",
        "Missing: ",
        paste(missing_vars, collapse = ", ")
      )
    }

    if (length(missing_vars) > 0 && !quiet) {
      message(
        "Note: Variables not in data (may be latent/derived): ",
        paste(missing_vars, collapse = ", ")
      )
    }

    # Handle id_col for matching to tree/structure
    row_ids <- NULL  # (re-initialised here for the flat path)
    if (is.data.frame(data)) {
      if (!is.null(id_col)) {
        if (!id_col %in% names(data)) {
          stop("id_col '", id_col, "' not found in data frame columns.")
        }
        row_ids <- data[[id_col]]
        # Remove id_col from variables to include (it's metadata, not a model variable)
        available_vars <- setdiff(available_vars, id_col)
      } else {
        # Try to use row names if they're meaningful (not just 1, 2, 3...)
        rn <- rownames(data)
        if (!is.null(rn) && !all(rn == as.character(seq_len(nrow(data))))) {
          row_ids <- rn
        }
      }
    }
    # Also check for variability-related columns (X_se, X_obs patterns)
    se_cols <- grep("_se$", names(data), value = TRUE)
    obs_cols <- grep("_obs$", names(data), value = TRUE)

    # Add generated dummy variables for categorical predictors
    dummy_vars <- character(0)
    if (!is.null(attr(data, "categorical_vars"))) {
      cat_vars <- attr(data, "categorical_vars")

      # Extract RHS variables (predictors) from fixed equations to filter dummies
      rhs_vars <- unique(unlist(lapply(equations, function(eq) {
        if (length(eq) == 3) {
          all.vars(eq[[3]])
        } else {
          character(0)
        }
      })))

      # Only generate dummies for variables used as predictors
      cat_vars <- cat_vars[names(cat_vars) %in% rhs_vars]

      dummy_vars <- unlist(lapply(cat_vars, function(x) x$dummies))
    }

    extra_cols <- c(se_cols, obs_cols, dummy_vars)

    # Convert to list format

    data_list <- list()
    for (var in c(available_vars, extra_cols)) {
      if (var %in% names(original_data)) {
        data_list[[var]] <- original_data[[var]]
      } else if (var %in% names(data)) {
        data_list[[var]] <- data[[var]]
      }
    }

    # Set names on vectors if we have row_ids
    if (!is.null(row_ids)) {
      for (var in names(data_list)) {
        if (
          is.vector(data_list[[var]]) &&
            length(data_list[[var]]) == length(row_ids)
        ) {
          names(data_list[[var]]) <- row_ids
        }
      }
    }

    # Preserve any existing attributes
    data_attrs <- attributes(original_data)

    # Combine data_list with random effect data
    data <- data_list
    if (length(random_data_updates) > 0) {
      for (nm in names(random_data_updates)) {
        if (!is.null(nm) && nm != "") {
          data[[nm]] <- random_data_updates[[nm]]
        }
      }
    }

    # Remove raw grouping variables from data passed to JAGS to avoid "Unused variable" warnings
    if (length(random_terms) > 0) {
      # Use setdiff to avoid error if variable not present (though they should be)
      vars_to_remove <- intersect(names(data), random_vars)
      if (length(vars_to_remove) > 0) {
        data[vars_to_remove] <- NULL
      }
    }

    if (!quiet) {
      message(
        "Converted data.frame to list with ",
        length(data),
        " variables: ",
        paste(names(data), collapse = ", ")
      )
    }
  } else {
    # If optimized hierarchical path was taken, we skipped the big block above.
    # We MUST merge random data updates (Prec matrices) now.
    if (is_hierarchical) {
      if (length(random_data_updates) > 0) {
        for (nm in names(random_data_updates)) {
          if (!is.null(nm) && nm != "") {
            data[[nm]] <- random_data_updates[[nm]]
          }
        }
      }
      # Restore attributes (categorical_vars) needed for model generation
      if (
        !is.null(hierarchical_info) &&
          !is.null(attr(hierarchical_info$data, "categorical_vars"))
      ) {
        attr(data, "categorical_vars") <- attr(
          hierarchical_info$data,
          "categorical_vars"
        )
      }
      data_attrs <- list()
    } else {
      # Standard Fallback
      # Ensure data is a list (crucial for adding matrices like VCV)
      # Preserve attributes (like categorical_vars) which are lost during as.list()
      data_attrs <- attributes(data)
      data <- as.list(data)
    }
  }

  # Compute polynomial values
  # We add them to original_data (for d-sep/residuals) but NOT to 'data' passed to JAGS
  # because JAGS creates deterministic nodes for them (var_pow2 <- var^2)
  if (!is.null(all_poly_terms)) {
    for (poly_term in all_poly_terms) {
      base_var <- poly_term$base_var
      power <- poly_term$power
      internal_name <- poly_term$internal_name

      # Check if base variable exists in data
      if (base_var %in% names(data)) {
        # Compute polynomial: x^2, x^3, etc.
        poly_vals <- data[[base_var]]^power

        # Add to original_data if possible
        if (is.list(original_data) || is.data.frame(original_data)) {
          original_data[[internal_name]] <- poly_vals
        }

        # If hierarchical, we MUST also add it to the source dataframes in hierarchical_info
        # Otherwise d-separation tests (which re-fetch data) won't find the new variable
        if (is_hierarchical && !is.null(hierarchical_info)) {
          # Find which level/dataframe holds the base variable
          for (lvl_name in names(hierarchical_info$data)) {
            if (base_var %in% names(hierarchical_info$data[[lvl_name]])) {
              # Add computed column to this level's dataframe
              source_df <- hierarchical_info$data[[lvl_name]]
              hierarchical_info$data[[lvl_name]][[internal_name]] <- source_df[[
                base_var
              ]]^power
              break
            }
          }
        }
      }
    }
  }

  # Restore categorical_vars if present
  if ("categorical_vars" %in% names(data_attrs)) {
    attr(data, "categorical_vars") <- data_attrs$categorical_vars
  }


  return(data)
}
