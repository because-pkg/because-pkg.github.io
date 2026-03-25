#' Format Data for Because Analysis
#'
#' Converts data from long format (one row per observation) to the list format
#' required by \code{\link{because}}.
#'
#' @param data A data.frame in long format with one row per observation.
#' @param species_col Name of the column containing species or unit identifiers (default: "SP").
#' @param tree A phylogenetic tree or other structural object. Optional. If provided, determines the order of units.
#'
#' @return A named list where each element is either:
#'   \itemize{
#'     \item A numeric vector (if all species have exactly 1 observation)
#'     \item A numeric matrix with species in rows and replicates in columns
#'   }
#'   Species are ordered to match \code{tree$tip.label} (if provided) or sorted alphabetically.
#'
#' @details
#' This function handles:
#' \itemize{
#'   \item Different numbers of replicates per species (creates rectangular matrix with NA padding)
#'   \item Missing values (NA)
#'   \item Automatic alignment with phylogenetic tree tip labels (if provided)
#' }
#'
#' When species have different numbers of replicates, the function creates a matrix
#' with dimensions (number of species) x (maximum number of replicates).
#' Species with fewer replicates are padded with NA values.
#'
#' If a tree is provided:
#' Species in the tree but not in the data will have all NA values.
#' Species in the data but not in the tree will be excluded with a warning.
#'
#' If no tree is provided:
#' All species in the data are included, sorted alphabetically by their ID.
#'
#' @examples
#' \dontrun{
#' # Example data in long format
#' data_long <- data.frame(
#'   SP = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp3"),
#'   BM = c(1.2, 1.3, 1.1, 2.1, 2.2, 1.8),
#'   NL = c(0.5, 0.6, NA, 0.7, 0.8, 0.9)
#' )
#'
#' # With tree
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   tree <- ape::read.tree(text = "(sp1:1,sp2:1,sp3:1);")
#'   data_list <- because_format_data(data_long, species_col = "SP", tree = tree)
#' }
#'
#' # Without tree (general repeated measures)
#' data_list_no_tree <- because_format_data(data_long, species_col = "SP")
#' }
#'
#' @export
because_format_data <- function(data, species_col = "SP", tree = NULL) {
    # Validate inputs
    if (!is.data.frame(data)) {
        stop("'data' must be a data.frame")
    }

    if (!species_col %in% names(data)) {
        stop(sprintf("Species column '%s' not found in data", species_col))
    }

    # No hard class check (allowed by get_order_labels_hook)

    # Get trait columns (everything except species column)
    trait_cols <- setdiff(names(data), species_col)

    if (length(trait_cols) == 0) {
        stop("No trait columns found in data")
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

    # Check for species alignment
    data_species <- unique(data[[species_col]])

    if (!is.null(tree)) {
        reference_labels <- get_order_labels_hook(tree)
        if (is.null(reference_labels)) {
            stop(
                "Could not extract labels from structure. Ensure it has tip labels or rownames."
            )
        }
        n_species <- length(reference_labels)

        missing_in_tree <- setdiff(data_species, reference_labels)
        missing_in_data <- setdiff(reference_labels, data_species)

        if (length(missing_in_tree) > 0) {
            warning(sprintf(
                "Species in data but not in tree (will be excluded): %s",
                paste(missing_in_tree, collapse = ", ")
            ))
            # Filter out species not in tree
            data <- data[data[[species_col]] %in% reference_labels, ]
        }

        if (length(missing_in_data) > 0) {
            message(sprintf(
                "Species in tree but not in data (will have NA values): %s",
                paste(missing_in_data, collapse = ", ")
            ))
        }
    } else {
        # No tree: use species in data, sorted
        reference_labels <- sort(data_species)
        n_species <- length(reference_labels)
    }

    # Count observations per species
    obs_counts <- table(data[[species_col]])
    max_reps <- max(obs_counts)

    # Initialize output list
    data_list <- list()

    # Process each trait
    for (trait in trait_cols) {
        # Create matrix: rows = species (in reference order), columns = replicates
        trait_matrix <- matrix(NA, nrow = n_species, ncol = max_reps)
        rownames(trait_matrix) <- reference_labels

        # Fill matrix with observations
        is_constant <- TRUE

        for (i in seq_along(reference_labels)) {
            sp <- reference_labels[i]

            if (sp %in% names(obs_counts)) {
                # Get observations for this species
                sp_values <- data[data[[species_col]] == sp, trait, drop = TRUE]
                n_obs <- length(sp_values)

                # Check for within-species variation
                # Stricter: if multiple values exist and they differ, it's not constant
                if (length(unique(na.omit(sp_values))) > 1) {
                    is_constant <- FALSE
                }

                # Fill in the observations (rest remain NA)
                trait_matrix[i, 1:n_obs] <- sp_values
            }
            # else: species not in data, row remains all NA
        }

        # If all species have <= 1 unique value (constant within species), convert to vector.
        # Usually, constant variables -> vector. Varying variables -> matrix.

        if (max_reps == 1 || is_constant) {
            # Take the first column (which contains the value for constant vars)
            # Note: for constant vars, all replicates are identical, so col 1 is sufficient.
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
