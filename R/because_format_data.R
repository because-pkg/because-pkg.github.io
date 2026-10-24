#' Format Long-Format Repeated Measures Data for Because Analysis
#'
#' Converts data from long format (multiple rows per unit/group, such as repeated
#' measurements per individual, site survey visits, or within-species replicates)
#' to the named list format (unit-level vectors and replicate matrices) used by
#' \code{\link{because}} for error-in-variables (\code{variability = "reps"}),
#' detection/occupancy models, and structured covariance models.
#'
#' @param data A \code{data.frame} in long format with one row per observation.
#' @param id_col Name of the column containing unit identifiers (e.g., individual,
#'   site, or species IDs; default: \code{"ID"}, falling back to \code{"SP"} if present).
#' @param structure Optional covariance structure object (e.g., a spatial/pedigree
#'   matrix with row names, a \code{because_structure} object, or a phylogenetic tree).
#'   If provided, determines the ordering of units in the output.
#' @param species_col Deprecated alias for \code{id_col} (retained for backward compatibility).
#' @param tree Deprecated alias for \code{structure} (retained for backward compatibility).
#'
#' @return A named list where each element is either:
#'   \itemize{
#'     \item A named numeric vector (if a variable is constant within each unit or all units have 1 observation)
#'     \item A numeric matrix with units in rows and replicates in columns (padded with \code{NA})
#'   }
#'   Units are ordered to match the labels of \code{structure} (if provided) or sorted alphabetically.
#'
#' @details
#' This function handles:
#' \itemize{
#'   \item Different numbers of replicates per unit (creates a rectangular matrix with \code{NA} padding)
#'   \item Automatic distinction between unit-constant covariates (returned as vectors) and replicate-varying measurements (returned as matrices)
#'   \item Missing values (\code{NA})
#'   \item Automatic alignment with covariance structure labels (e.g., matrix row names or phylogenetic tree tip labels, if \code{structure} is provided)
#' }
#'
#' When units have different numbers of replicates, the function creates a matrix
#' with dimensions \code{(number of units) x (maximum number of replicates)}.
#' Units with fewer replicates are padded with \code{NA} values.
#'
#' If a \code{structure} is provided:
#' \itemize{
#'   \item Units present in the structure but missing from \code{data} receive all \code{NA} values (enabling Bayesian imputation).
#'   \item Units present in \code{data} but absent from the structure are excluded with a warning.
#' }
#'
#' If no \code{structure} is provided:
#' All units in \code{data} are included, sorted alphabetically by their identifier.
#'
#' @examples
#' \dontrun{
#' # Example 1: General repeated measures / Error-in-Variables (no structure)
#' data_long <- data.frame(
#'   ID   = c("ind1", "ind1", "ind1", "ind2", "ind2", "ind3"),
#'   Age  = c(3, 3, 3, 5, 5, 2),          # Constant within individual -> vector
#'   Mass = c(1.2, 1.3, 1.1, 2.1, 2.2, 1.8), # Varying -> matrix (3 x 3)
#'   Cort = c(0.5, 0.6, NA, 0.7, 0.8, 0.9)   # Varying -> matrix (3 x 3)
#' )
#' data_list <- because_format_data(data_long, id_col = "ID")
#'
#' # Example 2: Aligning with a covariance matrix (e.g., spatial or pedigree)
#' prec_mat <- diag(3)
#' rownames(prec_mat) <- colnames(prec_mat) <- c("ind3", "ind1", "ind2")
#' data_aligned <- because_format_data(data_long, id_col = "ID", structure = prec_mat)
#'
#' # Example 3: Phylogenetic tree alignment (backward-compatible aliases)
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   tree <- ape::read.tree(text = "(ind1:1,ind2:1,ind3:1);")
#'   data_phylo <- because_format_data(data_long, species_col = "ID", tree = tree)
#' }
#' }
#'
#' @export
because_format_data <- function(
    data,
    id_col = "ID",
    structure = NULL,
    species_col = NULL,
    tree = NULL
) {
    # Validate inputs
    if (!is.data.frame(data)) {
        stop("'data' must be a data.frame")
    }

    # Backward compatibility: species_col -> id_col, tree -> structure
    if (!is.null(species_col)) {
        id_col <- species_col
    } else if (missing(id_col) && !("ID" %in% names(data)) && ("SP" %in% names(data))) {
        id_col <- "SP"
    }

    if (is.null(structure) && !is.null(tree)) {
        structure <- tree
    }

    if (is.null(id_col) || !id_col %in% names(data)) {
        stop(sprintf("Unit identifier column '%s' not found in data", as.character(id_col)))
    }

    # Get variable columns (everything except unit identifier column)
    trait_cols <- setdiff(names(data), id_col)

    if (length(trait_cols) == 0) {
        stop("No variable columns found in data")
    }

    # Detect and expand categorical variables (factor or character)
    categorical_vars <- list()
    for (col in trait_cols) {
        if (is.factor(data[[col]]) || is.character(data[[col]])) {
            # Get unique levels (excluding NA)
            if (is.factor(data[[col]])) {
                # Preserve factor ordering but exclude levels not present in data
                all_levels <- levels(data[[col]])
                present_levels <- unique(data[[col]][!is.na(data[[col]])])
                levels <- all_levels[all_levels %in% present_levels]
            } else {
                levels <- sort(unique(data[[col]][!is.na(data[[col]])]))
            }

            if (length(levels) < 2) {
                warning(sprintf(
                    "Categorical variable '%s' has fewer than 2 levels, skipping dummy coding",
                    col
                ))
                next
            }

            is_ord <- is.ordered(data[[col]])

            if (is_ord) {
                # Polynomial contrasts
                c_mat <- stats::contr.poly(length(levels))
                c_names <- colnames(c_mat)
                c_names[c_names == ".L"] <- "L"
                c_names[c_names == ".Q"] <- "Q"
                c_names[c_names == ".C"] <- "C"
                c_names <- gsub("\\^", "pow", c_names)

                dummies <- paste0(col, "_", c_names)

                categorical_vars[[col]] <- list(
                    levels = levels,
                    reference = "Polynomial Contrast",
                    dummies = dummies,
                    type = "ordered"
                )

                # Apply contrasts to data
                for (i in seq_along(dummies)) {
                    data[[dummies[i]]] <- c_mat[match(data[[col]], levels), i]
                }

                message(sprintf(
                    "Ordered categorical variable '%s' expanded to %d polynomial contrast(s)",
                    col,
                    length(levels) - 1
                ))
            } else {
                # Store categorical info (unordered treatment contrasts)
                categorical_vars[[col]] <- list(
                    levels = levels,
                    reference = levels[1],
                    dummies = paste0(col, "_", levels[-1]),
                    type = "unordered"
                )

                # Create dummy variables (K-1 for K levels)
                for (i in seq_along(levels)[-1]) {
                    dummy_name <- paste0(col, "_", levels[i])
                    data[[dummy_name]] <- as.numeric(data[[col]] == levels[i])
                }

                message(sprintf(
                    "Categorical variable '%s' expanded to %d dummy variable(s) | Reference: '%s'",
                    col,
                    length(levels) - 1,
                    levels[1]
                ))
            }
        }
    }

    # Update trait_cols: remove categoricals, add dummies
    if (length(categorical_vars) > 0) {
        original_categoricals <- names(categorical_vars)
        all_dummies <- unlist(lapply(categorical_vars, function(x) x$dummies))
        trait_cols <- c(setdiff(trait_cols, original_categoricals), all_dummies)
    }

    # Check for unit alignment
    data_units <- unique(as.character(data[[id_col]]))

    if (!is.null(structure)) {
        reference_labels <- get_order_labels_hook(structure)
        if (is.null(reference_labels)) {
            tree_obj <- get_tree_hook(structure)
            if (!is.null(tree_obj)) {
                reference_labels <- get_order_labels_hook(tree_obj)
            }
        }
        if (is.null(reference_labels)) {
            stop(
                "Could not extract labels from structure. Ensure it has tip labels or rownames."
            )
        }
        n_units <- length(reference_labels)

        missing_in_struct <- setdiff(data_units, reference_labels)
        missing_in_data <- setdiff(reference_labels, data_units)

        if (length(missing_in_struct) > 0) {
            warning(sprintf(
                "Units in data but not in structure (will be excluded): %s",
                paste(missing_in_struct, collapse = ", ")
            ))
            # Filter out units not in structure
            data <- data[data[[id_col]] %in% reference_labels, ]
        }

        if (length(missing_in_data) > 0) {
            message(sprintf(
                "Units in structure but not in data (will have NA values): %s",
                paste(missing_in_data, collapse = ", ")
            ))
        }
    } else {
        # No structure: use units in data, sorted alphabetically
        reference_labels <- sort(data_units)
        n_units <- length(reference_labels)
    }

    # Count observations per unit
    obs_counts <- table(as.character(data[[id_col]]))
    max_reps <- if (length(obs_counts) > 0) max(obs_counts) else 1L

    # Initialize output list
    data_list <- list()

    # Process each variable
    for (trait in trait_cols) {
        # Create matrix: rows = units (in reference order), columns = replicates
        trait_matrix <- matrix(NA, nrow = n_units, ncol = max_reps)
        rownames(trait_matrix) <- reference_labels

        # Fill matrix with observations
        is_constant <- TRUE

        for (i in seq_along(reference_labels)) {
            unit_id <- reference_labels[i]

            if (unit_id %in% names(obs_counts)) {
                # Get observations for this unit
                unit_values <- data[as.character(data[[id_col]]) == unit_id, trait, drop = TRUE]
                n_obs <- length(unit_values)

                # Check for within-unit variation
                # If multiple values exist and they differ, it's not constant
                if (length(unique(stats::na.omit(unit_values))) > 1) {
                    is_constant <- FALSE
                }

                # Fill in the observations (rest remain NA)
                trait_matrix[i, 1:n_obs] <- unit_values
            }
            # else: unit not in data, row remains all NA
        }

        # If all units have <= 1 unique value (constant within unit), convert to vector.
        # Constant variables -> vector. Varying variables -> matrix.
        if (max_reps == 1 || is_constant) {
            trait_vector <- trait_matrix[, 1]
            names(trait_vector) <- reference_labels
            data_list[[trait]] <- trait_vector
        } else {
            data_list[[trait]] <- trait_matrix
        }
    }
    # Store categorical variable info as attribute
    if (length(categorical_vars) > 0) {
        attr(data_list, "categorical_vars") <- categorical_vars
    }

    return(data_list)
}
