#' Auto-Detect Hierarchical Data Structure
#'
#' Infers levels, hierarchy, and link variables from a list of data.frames
#'
#' @param data List of data.frames at different hierarchical levels
#' @param eq_vars Character vector of variable names used in equations
#' @param quiet Logical, if TRUE suppress messages
#' @return List with 'levels', 'hierarchy', 'link_vars'
#' @keywords internal
auto_detect_hierarchical <- function(data, eq_vars, quiet = FALSE) {
    df_names <- names(data)

    # 1. Detect which equation variables are in which dataframe
    levels <- list()
    for (df_name in df_names) {
        # Skip if not a dataframe (e.g. user passed 'N=100' in data list)
        if (!is.data.frame(data[[df_name]])) {
            next
        }

        df <- data[[df_name]]
        vars_in_df <- intersect(eq_vars, colnames(df))
        if (length(vars_in_df) > 0) {
            levels[[df_name]] <- vars_in_df
        }
    }

    # Only trigger hierarchical mode if we have at least 2 levels (multiple dataframes)
    # Single dataframes in a list (common in JAGS/occupancy) are better handled
    # by standard list-to-dataframe conversion if they aren't truly hierarchical.
    if (length(levels) < 2) {
        if (!quiet && length(levels) == 1) {
            message(
                "Single level detected in list. Proceeding with standard joint data format."
            )
        }
        return(list(levels = NULL, hierarchy = NULL, link_vars = NULL))
    }

    # 2. Detect link variables (columns appearing in multiple dataframes)
    all_cols <- lapply(data, colnames)
    col_counts <- table(unlist(all_cols))
    link_vars <- names(col_counts[col_counts > 1])

    # Exclude equation variables from link_vars (they are data, not IDs)
    link_vars <- setdiff(link_vars, eq_vars)

    # 3. Infer hierarchy by row count (fewer rows = coarser level)
    row_counts <- sapply(data[names(levels)], nrow)
    ordered_levels <- names(sort(row_counts))
    hierarchy <- paste(ordered_levels, collapse = " > ")

    if (!quiet) {
        message("Auto-detected hierarchical structure:")
        for (lvl in ordered_levels) {
            vars <- levels[[lvl]]
            n <- nrow(data[[lvl]])
            message(sprintf(
                "  Level '%s' (N=%d): %s",
                lvl,
                n,
                paste(vars, collapse = ", ")
            ))
        }
        if (length(link_vars) > 0) {
            # Filter dummy variables (often end with _val or _spXX) to keep output clean
            # Heuristic: if var starts with another var in the list + "_", it's likely a dummy
            # Or simplified: if user asks to remove SpeciesID_spXX, we filter vars containing "_" if prefix exists?
            # Better: just look for standard dummy patterns if factor expansion happened.
            # But here we just want to suppress verbose output.

            display_links <- link_vars
            # Filter simple dummies (e.g. SpeciesID_sp1)
            is_dummy <- grepl("_.*[0-9]+$", display_links)
            # Only filter if we have non-dummies available
            if (any(!is_dummy)) {
                display_links <- display_links[!is_dummy]
            }

            message(sprintf(
                "  Link variables: %s",
                paste(display_links, collapse = ", ")
            ))
        }
        message(sprintf("  Hierarchy: %s", hierarchy))
    }

    return(list(
        levels = levels,
        hierarchy = hierarchy,
        link_vars = if (length(link_vars) > 0) link_vars else NULL
    ))
}


#' Validate Hierarchical Data Structure
#'
#' @param data List of data.frames at different hierarchical levels
#' @param levels List mapping variable names to level names
#' @param hierarchy Character string specifying nesting (e.g., "site_year > individual")
#' @param link_vars Character vector of variables that link levels
#' @keywords internal
validate_hierarchical_data <- function(
    data,
    levels,
    hierarchy,
    link_vars,
    latent_vars = NULL
) {
    # Check data is a named list of data.frames
    if (!is.list(data) || is.data.frame(data)) {
        stop(
            "For hierarchical data, 'data' must be a named list of data.frames"
        )
    }

    if (is.null(names(data)) || any(names(data) == "")) {
        stop("All elements of hierarchical 'data' must be named")
    }

    if (!all(sapply(data, is.data.frame))) {
        stop("All elements of hierarchical 'data' must be data.frames")
    }

    # Check levels is provided and valid
    if (is.null(levels)) {
        stop("'levels' argument required when using hierarchical data")
    }

    if (!is.list(levels) || is.null(names(levels))) {
        stop("'levels' must be a named list")
    }

    # Check that level names in 'levels' match dataset names in 'data'
    if (!all(names(levels) %in% names(data))) {
        missing <- setdiff(names(levels), names(data))
        stop("Level names not found in data: ", paste(missing, collapse = ", "))
    }

    # Check that all variables in levels exist in corresponding datasets
    # (Unless they are marked as latent)
    for (level_name in names(levels)) {
        vars <- levels[[level_name]]
        dataset <- data[[level_name]]

        # Filter out variables that are meant to be latent
        expected_vars <- setdiff(vars, latent_vars)
        missing_vars <- setdiff(expected_vars, colnames(dataset))

        if (length(missing_vars) > 0) {
            stop(
                "Variables not found in ",
                level_name,
                " dataset: ",
                paste(missing_vars, collapse = ", ")
            )
        }
    }

    # Check for variable overlap (variables should not appear in multiple levels)
    all_vars <- unlist(levels)
    if (any(duplicated(all_vars))) {
        dups <- all_vars[duplicated(all_vars)]
        stop(
            "Variables appear in multiple levels (not allowed): ",
            paste(unique(dups), collapse = ", ")
        )
    }

    # Validate link_vars if provided
    if (!is.null(link_vars)) {
        # Convert list to vector if needed
        link_vec <- if (is.list(link_vars)) unlist(link_vars) else link_vars

        # We cannot strictly enforce that EVERY link var is in EVERY dataset
        # because of multi-membership (e.g., Species and Site are independent levels)
        # Instead, verify that each dataset has AT LEAST one link variable,
        # UNLESS it's the very top level of a single hierarchy.
        # Actually, the safest validation is just ensuring they exist *somewhere* in the data
        all_cols <- unique(unlist(lapply(data, colnames)))
        missing_links <- setdiff(link_vec, all_cols)
        if (length(missing_links) > 0) {
            stop(
                "Link variables not found in any dataset: ",
                paste(missing_links, collapse = ", ")
            )
        }
    }

    # Validate hierarchy string if provided
    if (!is.null(hierarchy)) {
        if (!is.character(hierarchy) || length(hierarchy) != 1) {
            stop("'hierarchy' must be a single character string")
        }

        # Parse hierarchy (e.g., "site_year > individual" or "site > obs; species > obs")
        # First split by semicolon for multiple independent hierarchies
        hierarchy_paths <- strsplit(hierarchy, "\\s*;\\s*")[[1]]

        # Then extract all unique levels across all paths
        hierarchy_levels <- unique(unlist(lapply(hierarchy_paths, function(h) {
            strsplit(h, "\\s*>\\s*")[[1]]
        })))

        # Check all hierarchy levels exist in data
        missing_h <- setdiff(hierarchy_levels, names(data))
        if (length(missing_h) > 0) {
            stop(
                "Hierarchy levels not found in data: ",
                paste(missing_h, collapse = ", ")
            )
        }
    }

    invisible(TRUE)
}


#' Infer Variable Level
#'
#' Determine which hierarchical level a variable belongs to
#'
#' @param var Character, variable name
#' @param levels List mapping variable names to level names
#' @return Character, level name
#' @keywords internal
infer_variable_level <- function(var, levels) {
    for (level_name in names(levels)) {
        if (var %in% levels[[level_name]]) {
            return(level_name)
        }
    }

    # Variable not found in any level
    stop("Variable '", var, "' not found in any hierarchical level")
}


#' Get Data for Variables
#'
#' Determine the finest grain level needed for a set of variables
#' and return the appropriate dataset (with joining if needed)
#'
#' @param variables Character vector of variable names
#' @param data List of data.frames at different levels
#' @param levels List mapping variable names to level names
#' @param hierarchy Character string specifying nesting
#' @param link_vars Character vector of linking variables
#' @return data.frame
#' @keywords internal
get_data_for_variables <- function(
    variables,
    data,
    levels,
    hierarchy,
    link_vars
) {
    # If data is already a flat data.frame, just return the requested columns
    if (is.data.frame(data)) {
        # Return only the columns that exist in the flat data
        existing_vars <- intersect(variables, names(data))
        return(data[, existing_vars, drop = FALSE])
    }

    # Normalize link_vars: if passed as a named list (e.g. list(site="Site")), unlist to character vector
    if (is.list(link_vars)) {
        link_vars <- unlist(link_vars, use.names = FALSE)
    }

    # Determine which level each variable belongs to
    var_levels <- sapply(variables, function(v) {
        infer_variable_level(v, levels)
    })

    # Get unique levels needed
    needed_levels <- unique(var_levels)

    # If only one level, return that dataset
    if (length(needed_levels) == 1) {
        return(data[[needed_levels]])
    }

    # Multiple levels - need to determine which is finest grain and join
    # Parse hierarchy to get ordering (handling multi-membership like "site > obs; species > obs")
    hierarchy_paths <- strsplit(hierarchy, "\\s*;\\s*")[[1]]

    # Flatten to get a unified order for ranking depth (further right = finer grain)
    # We assign a depth score based on maximum index across all paths
    level_depths <- list()
    for (path in hierarchy_paths) {
        levels_in_path <- strsplit(path, "\\s*>\\s*")[[1]]
        for (i in seq_along(levels_in_path)) {
            lvl <- levels_in_path[i]
            if (is.null(level_depths[[lvl]]) || i > level_depths[[lvl]]) {
                level_depths[[lvl]] <- i
            }
        }
    }

    # Find the finest grain level among those needed (highest depth score)
    needed_depths <- sapply(needed_levels, function(l) level_depths[[l]] %||% 0)
    finest_idx <- which.max(needed_depths)
    finest_level <- needed_levels[finest_idx]

    # Start with finest grain dataset
    result <- data[[finest_level]]

    # Join data from coarser levels
    coarser_levels <- setdiff(needed_levels, finest_level)

    for (coarser_level in coarser_levels) {
        # Get variables from this level that we need
        vars_from_this_level <- names(var_levels)[var_levels == coarser_level]

        # Select those variables plus link vars from coarser dataset
        # Only include link_vars that actually exist in the coarser dataset
        valid_link_vars <- if (!is.null(link_vars)) {
            intersect(link_vars, names(data[[coarser_level]]))
        } else {
            character(0)
        }
        cols_to_select <- unique(c(valid_link_vars, vars_from_this_level))
        cols_to_select <- intersect(
            cols_to_select,
            names(data[[coarser_level]])
        )
        coarser_data <- data[[coarser_level]][,
            cols_to_select,
            drop = FALSE
        ]

        # Helper: Prune colliding variables from result that are claimed by coarser level
        # This ensures 'WINt' from enviro overrides 'WINt' in individual (if present but unmapped)
        # and prevents them from becoming unintended join keys.
        vars_to_add <- vars_from_this_level
        potential_collisions <- intersect(vars_to_add, names(result))

        # Don't prune if it's explicitly a link var (we want to use it for joining)
        if (!is.null(link_vars)) {
            potential_collisions <- setdiff(potential_collisions, link_vars)
        }

        if (length(potential_collisions) > 0) {
            # Drop colliding columns from result to prefer the coarser (source) version
            result[potential_collisions] <- NULL
        }

        # Join to result
        # The merge key must be columns that exist in BOTH result and coarser_data
        join_by <- intersect(names(result), names(coarser_data))
        # Optionally restrict to specified link_vars if provided
        if (!is.null(link_vars) && length(link_vars) > 0) {
            join_by_from_links <- intersect(link_vars, join_by)
            if (length(join_by_from_links) > 0) {
                join_by <- join_by_from_links
            }
        }

        if (length(join_by) == 0) {
            # No common columns - do a cross join (all combinations)
            result <- merge(result, coarser_data, all = FALSE)
        } else {
            result <- merge(result, coarser_data, by = join_by, all.x = TRUE)
        }
    }

    return(result)
}


#' Parse Hierarchy from Random Effects
#'
#' Extract hierarchical nesting structure from random effects formula
#'
#' @param random Formula specifying random effects
#' @param data List of data.frames at different hierarchical levels (optional, currently unused)
#' @return Character string hierarchy (e.g., "site_year > individual") or NULL
#' @keywords internal
parse_hierarchy_from_random <- function(random, data = NULL) {
    if (is.null(random)) {
        return(NULL)
    }

    # Convert to character
    random_str <- deparse(random)

    # Look for explicit nesting syntax (1|A/B)
    # After expansion, this will be (1|A) + (1|A:B)
    # We want to extract A > A:B pattern
    nesting_pattern <- "\\(1\\s*\\|\\s*([^/)]+)/([^)]+)\\)"

    if (grepl(nesting_pattern, random_str)) {
        matches <- regmatches(
            random_str,
            regexec(nesting_pattern, random_str)
        )[[1]]

        if (length(matches) >= 3) {
            # matches[2] = A, matches[3] = B
            # Hierarchy: A > A_B (Assuming B is nested in A)
            # Actually, (1|A/B) expands to (1|A) + (1|A:B).
            # So hierarchy is "A > B" (conceptually) or "A > id".
            # Let's just return "A > B" as a heuristic
            return(paste(matches[2], ">", matches[3]))
        }
    }

    return(NULL)
}

#' Auto-Stack Multispecies Data
#'
#' Converts a list of matrices (Wide) into a stacked dataframe (Long) for hierarchical analysis.
#' Replicates site and species covariates accordingly.
#'
#' @param data Input list of data (e.g. list(Y=list(Sp1=mat...), Hab=vec, Trait=vec))
#' @param equations List of model equations
#' @param quiet Logical, suppress messages
#'
#' @return List with components:
#'   \item{data}{New stacked dataframe}
#'   \item{random_part}{String to append to random formula (e.g. "+ (1|SpeciesID)")}
#'   \item{is_stacked}{Logical, whether stacking occurred}
#' @keywords internal
auto_stack_multispecies_data <- function(data, equations, quiet = FALSE) {
    if (!is.list(data) || is.data.frame(data)) {
        return(list(data = data, is_stacked = FALSE))
    }

    # 1. Identify "List Variables" that are in equations
    eq_vars <- unique(unlist(lapply(equations, all.vars)))
    list_vars <- names(data)[sapply(data, function(x) {
        is.list(x) && !is.data.frame(x)
    })]
    target_vars <- intersect(list_vars, eq_vars)

    if (length(target_vars) == 0) {
        return(list(data = data, is_stacked = FALSE))
    }

    # We assume the first target list variable drives the stacking (usually 'Y')
    # If multiple (e.g. Y and p_covs list), we assume they are aligned by name/order.
    pivot_var <- target_vars[1]
    pivot_list <- data[[pivot_var]]

    # Validate Pivot List
    if (length(pivot_list) < 2) {
        return(list(data = data, is_stacked = FALSE))
    } # Need at least 2 species to stack

    # Check if elements are matrices/dataframes
    valid_elems <- all(sapply(pivot_list, function(x) {
        is.matrix(x) || is.data.frame(x)
    }))
    if (!valid_elems) {
        return(list(data = data, is_stacked = FALSE))
    }

    # Get Dimensions
    # Assume all have same nrows (Sites)
    N_sites <- nrow(pivot_list[[1]])
    if (any(sapply(pivot_list, nrow) != N_sites)) {
        stop(
            "Auto-Stacking Error: All matrices in '",
            pivot_var,
            "' must have the same number of rows (Sites)."
        )
    }

    species_names <- names(pivot_list)
    if (is.null(species_names)) {
        species_names <- paste0("Sp", seq_along(pivot_list))
    }
    N_species <- length(species_names)

    if (!quiet) {
        message(
            "Auto-Stacking detected: Converting ",
            N_species,
            " species matrices into single hierarchical dataset."
        )
    }

    # --- Perform Stacking ---
    new_data <- list()

    # 1. Stack the Pivot Variable(s)
    # If Y is detections, we stack into a single large matrix (N_total x N_reps)
    # We need to ensure columns (reps) align too, or pad with NA?
    # For now, require same cols or use simplify2array logic.

    for (v in target_vars) {
        # Convert list of matrices to stacked matrix
        # Do simplistic stack: do.call(rbind, ...)
        # This stacks Sp1, then Sp2, etc. (Order: Species fast, or Site fast? Rbind is Species slow, Site fast)
        # Result: Rows 1..N = Sp1, Rows N+1..2N = Sp2.
        # This aligns with interaction: SpeciesID changes every N rows.

        mat_list <- lapply(data[[v]], as.matrix)
        # Check cols
        ref_cols <- ncol(mat_list[[1]])
        if (any(sapply(mat_list, ncol) != ref_cols)) {
            stop(
                "Auto-Stacking Error: All matrices must have same number of columns (Reps) for variable '",
                v,
                "'."
            )
        }
        stacked_mat <- do.call(rbind, mat_list)
        new_data[[v]] <- stacked_mat
    }

    # 2. Create ID Factors
    # SiteID: 1..N, 1..N, ...
    # SpeciesID: 1, 1.. (N times), 2, 2..

    new_data$SiteID <- factor(rep(1:N_sites, times = N_species))
    new_data$SpeciesID <- factor(
        rep(species_names, each = N_sites),
        levels = species_names
    )

    # 3. Handle Other Variables (Covariates)
    other_vars <- setdiff(names(data), target_vars)

    for (v in other_vars) {
        val <- data[[v]]

        # CASE A: Site Covariate (Vector of length N_sites)
        if (N_sites > 1 && length(val) == N_sites) {
            # Replicate for each species
            # Matches: 1..N (Sp1), 1..N (Sp2)...
            new_data[[v]] <- rep(val, times = N_species)
            if (!quiet) {
                message(
                    "  - Replicated site covariate '",
                    v,
                    "' (N=",
                    N_sites,
                    ")"
                )
            }
        } else if (N_species > 1 && length(val) == N_species) {
            # CASE B: Species Covariate (Vector of length N_species)
            # Replicate each element N_sites times
            # Matches: Sp1 (for 1..N), Sp2 (for 1..N)...

            # Check for name matching first
            if (!is.null(names(val)) && all(names(val) %in% species_names)) {
                # Match by name
                aligned_val <- val[species_names]
                new_data[[v]] <- rep(aligned_val, each = N_sites)
            } else {
                # Match by order
                new_data[[v]] <- rep(val, each = N_sites)
            }
            if (!quiet) {
                message(
                    "  - Replicated species trait '",
                    v,
                    "' (N=",
                    N_species,
                    ")"
                )
            }
        } else {
            # CASE C: Constant or Other?
            # If scalar, just keep as is? Or rep?
            # because usually handles scalars.
            # But dataframe construction requires vectors.
            # If length 1, rep to N_total
            if (length(val) == 1) {
                new_data[[v]] <- rep(val, N_sites * N_species)
            } else {
                warning(
                    "Variable '",
                    v,
                    "' has length ",
                    length(val),
                    " which does not match Sites (",
                    N_sites,
                    ") or Species (",
                    N_species,
                    "). Dropping from stacked data."
                )
            }
        }
    }

    # Convert to Data Frame
    # Note: Y might be a Matrix column. Data.frame supports matrix columns using I() or straight assignment if carefully done.
    # safely construct:
    # df <- data.frame(SiteID=..., SpeciesID=...)
    # df$Y <- stacked_mat

    df <- data.frame(SiteID = new_data$SiteID, SpeciesID = new_data$SpeciesID)

    for (v in names(new_data)) {
        if (v %in% c("SiteID", "SpeciesID")) {
            next
        }

        content <- new_data[[v]]
        if (is.matrix(content)) {
            # Assign matrix column
            df[[v]] <- content
        } else {
            df[[v]] <- content
        }
    }

    return(list(
        data = df,
        random_part = "+ (1|SpeciesID)",
        is_stacked = TRUE
    ))
}
