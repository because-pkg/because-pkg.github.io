#' Prepare Multiscale Data for JAGS
#'
#' Transforms multiscale dataframes into a flat list of vectors and ID indices suitable for JAGS
#'
#' @param hierarchical_info List containing 'data', 'levels', 'hierarchy', 'link_vars'
#' @param vars_needed Character vector of all variable names needed for the model
#' @return List with 'data_list' (for jags) and 'n_vec' (named vector of sample sizes)
#' @keywords internal
prepare_hierarchical_jags_data <- function(hierarchical_info, vars_needed) {
    data_list <- list()
    n_vec <- list()

    hierarchy_str <- hierarchical_info$hierarchy
    hierarchy_paths <- strsplit(hierarchy_str, "\\s*;\\s*")[[1]]

    unique_levels <- unique(unlist(lapply(hierarchy_paths, function(h) {
        trimws(strsplit(h, "\\s*>\\s*")[[1]])
    })))

    # 1. Extract variables from their native levels
    # Iterate through unique levels ensuring variables are extracted from their assigned level.
    for (lvl_name in unique_levels) {
        if (!lvl_name %in% names(hierarchical_info$data)) {
            next
        }

        df <- hierarchical_info$data[[lvl_name]]
        n_vec[[paste0("N_", lvl_name)]] <- nrow(df)

        # Identify variables that belong to this level
        lvl_vars <- hierarchical_info$levels[[lvl_name]]

        # Filter for only those needed in the model
        # Use known level variables + any other variables found in the dataframe (e.g. created dummies)
        vars_in_model <- unique(c(
            intersect(lvl_vars, vars_needed),
            intersect(names(df), vars_needed)
        ))

        for (v in vars_in_model) {
            # Check if variable is assigned to a specific level in hierarchy metadata
            # If so, only extract it from that level to avoid overwriting with wrong length (FK/PK issues)
            v_level <- NULL
            for (l_name in names(hierarchical_info$levels)) {
                if (v %in% hierarchical_info$levels[[l_name]]) {
                    v_level <- l_name
                    break
                }
            }

            if (!is.null(v_level) && v_level != lvl_name) {
                # Variable belongs to another level (e.g. Year belongs to Individual, don't extract from Enviro)
                next
            }

            if (!v %in% names(df)) {
                # Only error if it was EXPECTED here (i.e. this IS the assigned level)
                # If unassigned (e.g. sex_m), we expect it where found.
                if (!is.null(v_level) && v_level == lvl_name) {
                    stop(sprintf(
                        "Variable '%s' expected in level '%s' but not found in data.",
                        v,
                        lvl_name
                    ))
                }
                # If unassigned and missing here, just skip (look in other levels)
                next
            }
            data_list[[v]] <- df[[v]]
        }
    }

    # 2. Generate Indexing Vectors (Linking Child -> Parent)
    # Parse hierarchy (e.g., "Site > Species")
    # "Site" is child (finer), "Species" is parent (coarser) ??
    # Wait, usually hierarchy is written "Coarser > Finer" in my package?
    # Let's check: "Traits > Site_covs > Y"
    # User said: "Traits (50 rows) ... Site_covs (200 rows)"
    # So "Traits" is Top Level (Parent). "Site_covs" is Middle. "Y" is Bottom (Child).
    # Hierarchy string: "Traits > Site_covs > Y"
    # This implies Parent > Child.

    # 2. Generate Indexing Vectors (Linking Child -> Parent)
    for (path in hierarchy_paths) {
        levels_ordered <- trimws(strsplit(path, "\\s*>\\s*")[[1]])

        # Iterate strictly through adjacent pairs in the hierarchy
        # We need to link Level[k+1] (Child) to Level[k] (Parent)
        if (length(levels_ordered) > 1) {
            for (k in 1:(length(levels_ordered) - 1)) {
                parent_lvl <- levels_ordered[k]
                child_lvl <- levels_ordered[k + 1]

                parent_df <- hierarchical_info$data[[parent_lvl]]
                child_df <- hierarchical_info$data[[child_lvl]]

                # Determine link variable
                # We assume the link variable is the ID column of the PARENT
                # It must exist in both df
                common_cols <- intersect(names(parent_df), names(child_df))

                # Prefer 'link_vars' if specified
                if (!is.null(hierarchical_info$link_vars)) {
                    link_vec <- if (is.list(hierarchical_info$link_vars)) {
                        unlist(hierarchical_info$link_vars)
                    } else {
                        hierarchical_info$link_vars
                    }
                    link_col <- intersect(
                        common_cols,
                        link_vec
                    )
                    if (length(link_col) > 0) {
                        link_col <- link_col[1]
                    } else {
                        link_col <- NULL
                    }
                } else {
                    link_col <- if (length(common_cols) > 0) {
                        common_cols[1]
                    } else {
                        NULL
                    }
                }

                if (is.null(link_col)) {
                    stop(sprintf(
                        "Cannot link levels '%s' and '%s'. No common column found.",
                        parent_lvl,
                        child_lvl
                    ))
                }

                # Create Integer Index
                # Ensure factors match
                parent_vals <- as.character(parent_df[[link_col]])
                child_vals <- as.character(child_df[[link_col]])

                # Verify integrity
                missing_links <- setdiff(child_vals, parent_vals)
                if (length(missing_links) > 0) {
                    stop(sprintf(
                        "Level '%s' contains values in '%s' not found in parent level '%s'.",
                        child_lvl,
                        link_col,
                        parent_lvl
                    ))
                }

                idx_vec <- match(child_vals, parent_vals)

                # Naming convention: {Parent}_ID_in_{Child} ?? or just {Parent}_ID?
                # If we are in loop of Child, we access Parent[Parent_ID[i]]
                # But Child might not be the 'main' loop?
                # Actually, simple name: "{Parent}_ID" (vector of length Child)
                # Wait, if there are 3 levels A > B > C.
                # B needs A_ID (length N_B).
                # C needs B_ID (length N_C).
                # We name it specifically for the level it resides in?
                # No, JAGS models usually assume implicit context.
                # But here we have explicit loops.
                # Loop C uses: B_var[ B_index_in_C[i] ]
                # Loop B uses: A_var[ A_index_in_B[i] ]

                idx_name <- paste0(parent_lvl, "_idx_", child_lvl) # Explicit: Parent index vector residing in Child

                data_list[[idx_name]] <- idx_vec

                # Store mapping for model generator
                # We need to know: To access Parent from Child, use vector 'idx_name'
            }
        }

        # 3. Transitive Index Generation (Link Grandchild -> Grandparent)
        if (length(levels_ordered) > 2) {
            depth <- length(levels_ordered)
            for (gap in 2:(depth - 1)) {
                for (k in 1:(depth - gap)) {
                    parent_lvl <- levels_ordered[k]
                    child_lvl <- levels_ordered[k + gap]

                    # We need the intermediate level to bridge the gap
                    # We can use the level immediately preceding the child (k + gap - 1)
                    prev_child_lvl <- levels_ordered[k + gap - 1]

                    # We look for:
                    # 1. parent_idx_prev (Bridge -> Parent) - Created in previous gap iteration or adjacent step
                    # 2. prev_idx_child (Child -> Bridge)   - Created in adjacent step

                    parent_vec_name <- paste0(
                        parent_lvl,
                        "_idx_",
                        prev_child_lvl
                    )
                    prev_vec_name <- paste0(prev_child_lvl, "_idx_", child_lvl)

                    if (is.null(data_list[[parent_vec_name]])) {
                        # Should not occur if hierarchy is correct
                        warning(sprintf(
                            "Missing bridge index '%s' for transitive link.",
                            parent_vec_name
                        ))
                        next
                    }
                    if (is.null(data_list[[prev_vec_name]])) {
                        warning(sprintf(
                            "Missing child index '%s' for transitive link.",
                            prev_vec_name
                        ))
                        next
                    }

                    parent_vec <- data_list[[parent_vec_name]]
                    prev_vec <- data_list[[prev_vec_name]]

                    # Chain them: Parent[i] = Parent[ Prev[ Child[i] ] ]
                    # Effectively: parent_vec[ prev_vec ]

                    transitive_vec <- parent_vec[prev_vec]

                    transitive_name <- paste0(parent_lvl, "_idx_", child_lvl)
                    data_list[[transitive_name]] <- transitive_vec
                }
            }
        }
    }

    # 4. Generate Level-Specific Zeros Vectors (for structured residuals)
    for (lvl_name in unique_levels) {
        n_val <- n_vec[[paste0("N_", lvl_name)]]
        if (!is.null(n_val)) {
            data_list[[paste0("zeros_", lvl_name)]] <- rep(0, n_val)
        }
    }

    return(list(data_list = data_list, n_vec = n_vec))
}
