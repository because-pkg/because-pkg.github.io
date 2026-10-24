#' Auto-Detect Multiscale Data Structure
#'
#' Infers scales, hierarchy, and link variables from a list of data.frames
#'
#' @param data List of data.frames at different scales
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
        message("Auto-detected multiscale structure:")
        for (lvl in ordered_levels) {
            vars <- levels[[lvl]]
            n <- nrow(data[[lvl]])
            message(sprintf(
                "  Scale '%s' (N=%d): %s",
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


#' Auto-Detect Structure Levels
#'
#' Maps covariance structures (e.g. phylo, spatial) to hierarchical levels
#' by matching dimensions and ID values.
#'
#' @param structure List of structure objects
#' @param hierarchical_info List containing 'data', 'levels', etc.
#' @param quiet Logical
#' @return Named list mapping structure names to level names
#' @keywords internal
auto_detect_structure_levels <- function(structure, hierarchical_info, quiet = FALSE) {
    if (is.null(structure) || is.null(hierarchical_info)) {
        return(NULL)
    }

    levels <- hierarchical_info$levels
    data_list <- hierarchical_info$data
    s_levels <- list()

    for (s_name in names(structure)) {
        s_obj <- structure[[s_name]]

        # 1. Try to detect dimension of structure
        n_struct <- 0
        names_struct <- NULL

        if (inherits(s_obj, "phylo")) {
            n_struct <- length(s_obj$tip.label)
            names_struct <- s_obj$tip.label
        } else if (is.matrix(s_obj)) {
            n_struct <- nrow(s_obj)
            names_struct <- rownames(s_obj)
        } else if (inherits(s_obj, "multiPhylo")) {
            n_struct <- length(s_obj[[1]]$tip.label)
            names_struct <- s_obj[[1]]$tip.label
        }

        if (n_struct == 0) {
            next
        }
        
        # 2. Match against levels
        best_lvl <- NULL
        for (lvl_name in names(data_list)) {
            df <- data_list[[lvl_name]]
            if (!is.data.frame(df)) next
            n_lvl <- nrow(df)
            
            # Match by sample size
            if (n_lvl == n_struct) {
                # Match by names if possible
                if (!is.null(names_struct)) {
                    # Find potential ID columns in this level
                    id_cols <- names(df)[sapply(df, function(x) !is.numeric(x) || all(x %% 1 == 0, na.rm = TRUE))]
                    
                    match_found <- FALSE
                    cat_meta <- attr(data_list, "categorical_vars")
                    
                    for (col in id_cols) {
                        df_vals <- as.character(df[[col]])
                        
                        # [CATEGORICAL MAPPING] If encoded as integers, map back to labels
                        if (!is.null(cat_meta) && col %in% names(cat_meta)) {
                            lvls <- cat_meta[[col]]$levels
                            df_vals <- lvls[as.integer(df_vals)]
                        }

                        if (all(names_struct %in% df_vals)) {
                            match_found <- TRUE
                            break
                        }
                    }
                    
                    if (match_found) {
                        best_lvl <- lvl_name
                        break
                    }
                } else {
                    # Fallback to just sample size if no names provided in structure (e.g. anonymous distance matrix)
                    # We only accept this if it's the ONLY level with this sample size to avoid ambiguity.
                    all_n <- sapply(data_list, nrow)
                    if (sum(all_n == n_struct) == 1) {
                        best_lvl <- lvl_name
                        break
                    }
                }
            }
        }

        if (!is.null(best_lvl)) {
            s_levels[[s_name]] <- best_lvl
            if (!quiet) {
                message(sprintf("Auto-detected structure level: '%s' -> '%s'", s_name, best_lvl))
            }
        }
    }

    if (length(s_levels) == 0) {
        return(NULL)
    }
    return(s_levels)
}


#' Validate Multiscale Data Structure
#'
#' @param data List of data.frames at different scales
#' @param levels List mapping variable names to level names
#' @param hierarchy Character string specifying nesting (e.g., "site_year > individual")
#' @param link_vars Character vector of variables that link levels
#' @keywords internal
validate_hierarchical_data <- function(
    data,
    levels,
    hierarchy,
    link_vars,
    latent_vars = NULL,
    equations = NULL,
    deterministic_vars = NULL
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
        cat("\n[DEBUG] validate_hierarchical_data failing!\n")
        cat("  names(levels): ", paste(names(levels), collapse = ", "), "\n")
        cat("  names(data): ", paste(names(data), collapse = ", "), "\n")
        cat("  missing: ", paste(missing, collapse = ", "), "\n\n")
        stop("Level names not found in data: ", paste(missing, collapse = ", "))
    }

    # Combine provided latent variables with deterministic variables
    # for the validation skip list.
    latent_vars <- unique(c(latent_vars, deterministic_vars))

    # Check that all variables in levels exist in corresponding datasets
    # (Unless they are marked as latent or deterministic)
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
#' Determine which multiscale level a variable belongs to
#'
#' @param equations Optional list of model formulas for context
#' @param latent Character vector of latent variables
#' @param hierarchy Optional hierarchy string (e.g. "individual > obs") to determine ordering
#' @return Character, level name
#' @keywords internal
infer_variable_level <- function(var, levels, data = NULL, equations = NULL, latent = NULL, hierarchy = NULL) {
    # 1. Check explicit levels mapping
    for (level_name in names(levels)) {
        if (var %in% levels[[level_name]]) {
            return(level_name)
        }
    }

    # 2. Search data columns if not in levels list (fallback for structural IDs)
    if (!is.null(data)) {
        for (level_name in names(data)) {
            if (is.data.frame(data[[level_name]]) && var %in% colnames(data[[level_name]])) {
                return(level_name)
            }
        }
    }

    # 3. [AUTO-INFERENCE] If variable is latent, infer level from predictors
    if (!is.null(latent) && var %in% latent && !is.null(equations)) {
        # Find which equation this variable is the response of
        for (eq in equations) {
            resp <- as.character(formula(eq))[2]
            if (resp == var) {
                # Success! Infer level from predictors
                preds <- all.vars(formula(eq))[-1]
                if (length(preds) > 0) {
                    # Get levels of all predictors
                    pred_lvls <- na.omit(vapply(preds, function(p) {
                        # Recursive call WITHOUT passing equations/latent to avoid infinite loop
                        # (We only infer one level deep for simplicity)
                        tryCatch(infer_variable_level(p, levels, data), error = function(e) NA_character_)
                    }, character(1)))
                    
                    if (length(pred_lvls) > 0) {
                        if (!is.null(hierarchy)) {
                            # Rank levels according to hierarchy (coarsest first)
                            h_lvls <- trimws(strsplit(hierarchy, "\\s*>\\s*")[[1]])
                            found <- intersect(h_lvls, pred_lvls)
                            if (length(found) > 0) return(found[1])
                        }
                        # Fallback to the first found level
                        return(pred_lvls[1])
                    }
                }
            }
        }
        
        # If it's a latent parent (not a response), find its children
        for (eq in equations) {
            preds <- all.vars(formula(eq))[-1]
            if (var %in% preds) {
                # Infer level from the response variable of this equation
                resp <- as.character(formula(eq))[2]
                resp_lvl <- tryCatch(infer_variable_level(resp, levels, data), error = function(e) NA_character_)
                if (!is.na(resp_lvl)) return(resp_lvl)
            }
        }
    }

    # Variable not found in any level
    stop("Variable '", var, "' not found in any hierarchical level")
}

#' Parse Hierarchy from Random Effects
#'
#' Extract multiscale nesting structure from random effects formula
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

#' Get Depth of a Level in Hierarchy
#'
#' Correctly handles multi-chain hierarchies (semicolons).
#' Returns the 1-indexed depth (1 = coarsest).
#' If a level appears in multiple chains, returns the maximum depth.
#'
#' @param lvl Level name to find
#' @param hierarchy_str Hierarchy string
#' @return Integer depth
#' @keywords internal
get_level_depth <- function(lvl, hierarchy_str) {
  if (is.null(lvl) || is.null(hierarchy_str)) return(NA)
  
  paths <- strsplit(hierarchy_str, "\\s*;\\s*")[[1]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1]]))
  
  depths <- sapply(paths, function(p) {
    idx <- match(lvl, p)
    if (is.na(idx)) return(NA)
    return(idx)
  })
  
  if (all(is.na(depths))) return(NA)
  return(max(depths, na.rm = TRUE))
}
