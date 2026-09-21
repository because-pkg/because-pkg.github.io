#' Preprocess data and arguments for because()
#'
#' @param data The input data frame or list
#' @param equations The equations list
#' @param distribution The deprecated distribution arg
#' @param family The family arg
#' @param structure The structure arg
#' @return A list containing cleaned variables
#' @noRd
preprocess_because_args <- function(data, equations, distribution, family, structure) {
  # Input validation
  if (is.null(data)) {
    stop("Argument 'data' must be provided.")
  }
  if (is.null(equations)) {
    stop("Argument 'equations' must be provided.")
  }

  # --- Clean Data: Ensure matrix (e.g. from scale()) is a data.frame ---
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # --- Clean Data: Coerce 1D matrices (e.g., from scale()) to vectors ---
  scale_info <- list()
  if (is.data.frame(data)) {
    for (nm in names(data)) {
      if (!is.null(attr(data[[nm]], "scaled:center"))) {
        scale_info[[nm]] <- list(
          center = attr(data[[nm]], "scaled:center"),
          scale = attr(data[[nm]], "scaled:scale")
        )
      }
      if (is.matrix(data[[nm]]) && ncol(data[[nm]]) == 1) {
        data[[nm]] <- as.numeric(data[[nm]])
      }
    }
  } else if (is.list(data)) {
    for (df_name in names(data)) {
      if (is.data.frame(data[[df_name]])) {
        for (nm in names(data[[df_name]])) {
          if (!is.null(attr(data[[df_name]][[nm]], "scaled:center"))) {
            scale_info[[nm]] <- list(
              center = attr(data[[df_name]][[nm]], "scaled:center"),
              scale = attr(data[[df_name]][[nm]], "scaled:scale")
            )
          }
          if (is.matrix(data[[df_name]][[nm]]) && ncol(data[[df_name]][[nm]]) == 1) {
            data[[df_name]][[nm]] <- as.numeric(data[[df_name]][[nm]])
          }
        }
      }
    }
  }

  # --- Backward Compatibility: distribution -> family ---
  if (!is.null(distribution)) {
    warning(
      "Argument 'distribution' is deprecated and will be removed in future versions. Please use 'family' instead."
    )
    if (is.null(family)) {
      family <- distribution
    }
  }

  # --- Family Object Normalization (Custom Families Support) ---
  family_objects <- list()
  if (!is.null(family)) {
    if (inherits(family, "because_family")) {
      family_objects[["_default"]] <- family
      family <- setNames(family$family, "_default")
    } else if (is.list(family) && !is.null(names(family))) {
      for (nm in names(family)) {
        if (inherits(family[[nm]], "because_family")) {
          family_objects[[nm]] <- family[[nm]]
          family[[nm]] <- family[[nm]]$family
        }
      }
      family <- unlist(family)
    }
    
    supported_families <- c(
      "gaussian", "normal", "poisson", "negbinomial", "zip", "zinb", 
      "binomial", "bernoulli", "multinomial", "ordinal"
    )
    
    for (node_name in names(family)) {
      specified_dist <- tolower(family[[node_name]])
      if (!(specified_dist %in% supported_families)) {
        stop(
          sprintf(
            "Error: The distribution family '%s' specified for node '%s' is not supported in 'because'. Supported families are: %s.", 
            specified_dist, node_name, paste(supported_families, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  }

  family_obj <- family
  if (!is.null(family) && ("occupancy" %in% family || "cmr" %in% family)) {
    class(family_obj) <- c("because_family_occupancy", class(family_obj))
  }
  structure_obj <- structure
  if (!is.null(structure)) {
    class(structure_obj) <- c("because_structure", class(structure_obj))
  }

  equations <- normalize_equations_hook(family_obj, equations, data = data)

  return(list(
    data = data,
    scale_info = scale_info,
    family = family,
    family_objects = family_objects,
    family_obj = family_obj,
    structure_obj = structure_obj,
    equations = equations
  ))
}

#' Preprocess categorical variables (character/factor) to integer codes and dummies
#'
#' @param data A data.frame or list of data.frames
#' @param quiet Logical; whether to suppress informational messages
#' @return The modified data object with categorical_vars attribute
#' @keywords internal
preprocess_categorical_vars <- function(
  data,
  target_vars = NULL,
  dummy_vars = NULL,
  exclude_cols = NULL,
  quiet = FALSE,
  expand_ordered = FALSE
) {
  if (is.null(data)) {
    return(NULL)
  }

  # Recursively process lists of data frames (hierarchical data)
  if (is.list(data) && !is.data.frame(data)) {
    all_cat_vars <- list()

    # Process each element (usually levels in hierarchy)
    for (i in seq_along(data)) {
      if (is.data.frame(data[[i]])) {
        processed <- preprocess_categorical_vars(
          data[[i]],
          target_vars = target_vars,
          dummy_vars = dummy_vars,
          exclude_cols = exclude_cols,
          quiet = TRUE, # Suppress output for recursive calls to avoid noise
          expand_ordered = expand_ordered
        )
        data[[i]] <- processed
        # Collect categorical vars metadata
        level_cat_vars <- attr(processed, "categorical_vars")
        if (!is.null(level_cat_vars)) {
          # Use utils::modifyList if available, or manual merge
          for (name in names(level_cat_vars)) {
            all_cat_vars[[name]] <- level_cat_vars[[name]]
          }
        }
      }
    }

    # Attach merged metadata to the top-level list
    attr(data, "categorical_vars") <- all_cat_vars
    return(data)
  }

  # Process single data frame
  if (!is.data.frame(data)) {
    return(data)
  }

  # Determine which columns to check
  check_cols <- if (!is.null(target_vars)) {
    intersect(names(data), target_vars)
  } else {
    names(data)
  }

  if (length(check_cols) == 0) {
    return(data)
  }

  # Detect categorical variables: either by type or if they already have metadata
  cat_metadata <- attr(data, "categorical_vars")
  char_cols <- sapply(check_cols, function(col) {
    x <- data[[col]]
    is.character(x) || is.factor(x) || (!is.null(cat_metadata) && col %in% names(cat_metadata))
  })

  if (any(char_cols)) {
    categorical_vars <- list()
    if (!is.null(attr(data, "categorical_vars"))) {
      categorical_vars <- attr(data, "categorical_vars")
    }

    col_names <- check_cols[char_cols]
    for (col in col_names) {
      if (!is.null(exclude_cols) && col %in% exclude_cols) {
        next
      }

      # [FIX] Preserve existing categorical metadata if present (to keep levels consistent in subsets)
      existing_metadata <- categorical_vars[[col]]
      if (!is.null(existing_metadata)) {
        levels <- existing_metadata$levels
        is_ord_prev <- !is.null(existing_metadata$type) && existing_metadata$type == "ordered"
        
        # If already integer/numeric, it's likely already encoded from a previous call.
        # Use integer-to-label mapping to avoid wiping the data with NAs.
        if (is.numeric(data[[col]])) {
          f_vals <- factor(data[[col]], levels = seq_along(levels), labels = levels, ordered = is_ord_prev)
        } else {
          # Use existing levels to ensure integer encoding (1, 2, 3...) remains consistent
          # during d-separation tests on data subsets.
          f_vals <- factor(data[[col]], levels = levels, ordered = is_ord_prev)
        }
      } else {
        # New variable: Convert to factor to discover levels
        f_vals <- factor(data[[col]])
        levels <- levels(f_vals)
      }

      if (length(levels) < 2) {
        if (!quiet) {
          warning(sprintf(
            "Variable '%s' has < 2 levels. Converting to numeric constant.",
            col
          ))
        }
        data[[col]] <- as.numeric(f_vals)
      } else {
        # Store metadata for model expansion
        is_ord <- is.ordered(f_vals)

        if (is.null(existing_metadata)) {
          if (is_ord) {
            if (expand_ordered) {
              c_mat <- stats::contr.poly(length(levels))
              if (ncol(c_mat) > 2) {
                c_mat <- c_mat[, 1:2, drop = FALSE]
              }
              c_names <- colnames(c_mat)
              c_names[c_names == ".L"] <- "L"
              c_names[c_names == ".Q"] <- "Q"
              c_names[c_names == ".C"] <- "C"
              c_names <- gsub("\\^", "pow", c_names)

              categorical_vars[[col]] <- list(
                levels = levels,
                reference = "Polynomial Contrast",
                dummies = paste0(col, "_", c_names),
                type = "ordered",
                contrasts = c_mat
              )
            } else {
              # Default: Linear only (centered integers)
              categorical_vars[[col]] <- list(
                levels = levels,
                reference = "Numeric (Centered)",
                dummies = paste0(col, "_L"),
                type = "ordered"
              )
            }
          } else {
            categorical_vars[[col]] <- list(
              levels = levels,
              reference = levels[1],
              dummies = paste0(col, "_", levels[-1]),
              type = "unordered"
            )
          }
        } else {
          # Update existing metadata if needed
          if (is_ord && is.null(existing_metadata$contrasts) && expand_ordered) {
            c_mat <- stats::contr.poly(length(levels))
            if (ncol(c_mat) > 2) {
              c_mat <- c_mat[, 1:2, drop = FALSE]
            }
            c_names <- colnames(c_mat)
            c_names[c_names == ".L"] <- "L"
            c_names[c_names == ".Q"] <- "Q"
            c_names[c_names == ".C"] <- "C"
            c_names <- gsub("\\^", "pow", c_names)

            categorical_vars[[col]]$reference <- "Polynomial Contrast"
            categorical_vars[[col]]$dummies <- paste0(col, "_", c_names)
            categorical_vars[[col]]$type <- "ordered"
            categorical_vars[[col]]$contrasts <- c_mat
          } else if (is_ord && is.null(existing_metadata$dummies) && !expand_ordered) {
            categorical_vars[[col]]$reference <- "Numeric (Centered)"
            categorical_vars[[col]]$dummies <- paste0(col, "_L")
            categorical_vars[[col]]$type <- "ordered"
          }
        }

        # Convert to integer codes for JAGS
        data[[col]] <- as.integer(f_vals)

        # Generate Dummy Variables explicitly ONLY if requested
        if (is.null(dummy_vars) || col %in% dummy_vars) {
          dummies <- categorical_vars[[col]]$dummies

          if (!quiet && length(levels) > 500) {
            message(sprintf(
              "Generating %d dummy variables for '%s'... this may take a moment.",
              length(levels) - 1,
              col
            ))
          }
          if (is_ord) {
            if (expand_ordered) {
              c_mat <- categorical_vars[[col]]$contrasts
              for (k in seq_along(dummies)) {
                data[[dummies[k]]] <- c_mat[data[[col]], k]
              }
            } else {
              # Linear centered (dummy name is col_L)
              raw_codes <- data[[col]]
              data[[paste0(col, "_L")]] <- raw_codes - mean(raw_codes, na.rm = TRUE)
            }
          } else {
            for (k in 2:length(levels)) {
              lev_name <- levels[k]
              dummy_col_name <- dummies[k - 1]
              # Create binary column: 1 if matches this level, 0 if not (NAs remain NA)
              data[[dummy_col_name]] <- as.numeric(data[[col]] == k)
            }
          }
        }

        if (!quiet) {
          message(sprintf(
            "Converted categorical '%s' to integers (1..%d). Reference: '%s'%s",
            col,
            length(levels),
            levels[1],
            if (is.null(dummy_vars) || col %in% dummy_vars) {
              ""
            } else {
              " (Dummies skipped)"
            }
          ))
        }
      }
    }
    attr(data, "categorical_vars") <- categorical_vars
  }

  return(data)
}
