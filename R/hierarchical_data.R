
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
    link_vars,
    equations = NULL,
    latent = NULL
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
        infer_variable_level(v, levels, data = data, equations = equations, latent = latent, hierarchy = hierarchy)
    })

    # Get unique levels needed
    needed_levels <- unique(var_levels)

    # If only one level, return that dataset (Correct Resolution!)
    if (length(needed_levels) == 1) {
        return(data[[needed_levels]])
    }

    # Multiple levels - need to determine which is locally finest grain and join
    # Parse hierarchy to get ordering (handling multi-membership like "site > obs; species > obs")
    hierarchy_paths <- strsplit(hierarchy, "\\s*;\\s*")[[1]]

    # Flatten to get a unified order for ranking depth (further right = finer grain)
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

    # Find the finest grain level among those *actually* needed for this sub-test
    needed_depths <- sapply(needed_levels, function(l) level_depths[[l]] %||% 0)
    finest_idx <- which.max(needed_depths)
    finest_level_for_test <- needed_levels[finest_idx]

    # Start with THIS test's finest grain dataset (Resolution Locking!)
    result <- data[[finest_level_for_test]]

    # Join data from coarser levels
    coarser_levels <- setdiff(needed_levels, finest_level_for_test)

    # Sort coarser levels by depth (Descending: Join immediate parents FIRST to ensure keys are available)
    # e.g., for site > survey > obs, join survey to obs first, then site.
    coarser_depths <- sapply(coarser_levels, function(l) level_depths[[l]] %||% 0)
    coarser_levels <- coarser_levels[order(coarser_depths, decreasing = TRUE)]

    for (coarser_level in coarser_levels) {
        # Get variables from this level that we need
        vars_from_this_level <- names(var_levels)[var_levels == coarser_level]

        # Select those variables plus link vars from coarser dataset
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

        # Prune colliding variables (prefer coarser source)
        vars_to_add <- vars_from_this_level
        potential_collisions <- intersect(vars_to_add, names(result))
        if (!is.null(link_vars)) {
            potential_collisions <- setdiff(potential_collisions, link_vars)
        }
        if (length(potential_collisions) > 0) {
            result[potential_collisions] <- NULL
        }

        # Safe Keyed Join
        join_by <- intersect(names(result), names(coarser_data))
        if (!is.null(link_vars) && length(link_vars) > 0) {
            join_by_from_links <- intersect(link_vars, join_by)
            if (length(join_by_from_links) > 0) {
                join_by <- join_by_from_links
            }
        }

        if (length(join_by) == 0) {
            # DANGEROUS: Cross-Hierarchy Join (e.g. Species vs Site)
            # This happens when two variables are in different branches of the hierarchy.
            # We must warn and perform the merge but minimize the Cartesian impact by checking
            # if we can find any ID column in data[[lvl]].
            # For your manuscript: Abundance ~ Site will work because Abundance has SiteID.
            # But Species ~ Site is orthogonal.
            result <- merge(result, coarser_data, all = FALSE)
        } else {
            result <- merge(result, coarser_data, by = join_by, all.x = TRUE)
        }
    }

    return(result)
}


#' Auto-Stack Multispecies Data
#'
#' Converts a list of matrices (Wide) into a stacked dataframe (Long) for multiscale analysis.
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
#' Aggregate Multiscale Data to a Target Resolution
#'
#' Collapses finer-level variables to a coarser grain by taking the mean.
#' Used for cross-scale d-separation tests where tests must be locked to the
#' resolution of the coarsest variable involved.
#'
#' @param data List of dataframes
#' @param levels List mapping variables to levels
#' @param hierarchy Hierarchy string (e.g. "Year > Individual")
#' @param target_lvl The resolution to aggregate everything to
#' @return List of dataframes truncated/aggregated to stop at target_lvl
#' @keywords internal
aggregate_multiscale_data <- function(data, levels, hierarchy, target_lvl) {
  if (is.null(data) || is.null(target_lvl)) return(data)
  if (is.data.frame(data)) return(data) # Handle flat case (no-op or error?)

  # 1. Identify levels to keep and levels to aggregate
  all_lvls <- unique(unlist(lapply(strsplit(hierarchy, "\\s*[>;]\\s*"), trimws)))
  t_depth <- get_level_depth(target_lvl, hierarchy)
  if (is.na(t_depth)) return(data)

  keep_lvls <- Filter(function(l) {
      d <- get_level_depth(l, hierarchy)
      !is.na(d) && d <= t_depth
  }, all_lvls)
  agg_lvls <- Filter(function(l) {
      d <- get_level_depth(l, hierarchy)
      !is.na(d) && d > t_depth
  }, all_lvls)

  if (length(agg_lvls) == 0) return(data[names(data) %in% keep_lvls])

  # 2. Start with keep levels
  new_data <- data[names(data) %in% keep_lvls]
  target_df <- new_data[[target_lvl]]

  # 3. Aggregate each finer level to the target level
  for (lvl in agg_lvls) {
    vars_to_agg <- levels[[lvl]]
    if (length(vars_to_agg) == 0) next

    # Find link column (ID of target_lvl in lvl dataframe)
    link_col <- target_lvl
    if (!(link_col %in% names(data[[lvl]]))) {
        # Fallback to shared columns
        common <- intersect(names(data[[lvl]]), names(target_df))
        if (length(common) > 0) link_col <- common[1]
    }

    if (!(link_col %in% names(data[[lvl]]))) next

    # Only aggregate numeric variables
    is_numeric <- vapply(data[[lvl]][, vars_to_agg, drop=FALSE], is.numeric, logical(1))
    numeric_vars <- vars_to_agg[is_numeric]

    if (length(numeric_vars) > 0) {
        agg_df <- stats::aggregate(
          data[[lvl]][, numeric_vars, drop = FALSE],
          by = list(ID = data[[lvl]][[link_col]]),
          FUN = mean,
          na.rm = TRUE
        )
        names(agg_df)[1] <- link_col

        # Merge into target_df
        target_df <- merge(target_df, agg_df, by = link_col, all.x = TRUE)
    }
  }

  new_data[[target_lvl]] <- target_df
  return(new_data)
}

flatten_for_python <- function(data_list, link_vars = NULL) {
    if (is.data.frame(data_list) || (is.list(data_list) && !any(sapply(data_list, is.data.frame)))) {
        flat_data <- as.list(data_list)
    } else {
        flat_data <- list()
        for (lvl_name in names(data_list)) {
            if (is.data.frame(data_list[[lvl_name]])) {
                for (col in names(data_list[[lvl_name]])) {
                    flat_data[[col]] <- data_list[[lvl_name]][[col]]
                }
            }
        }
        
        if (!is.null(link_vars)) {
            for (link in unlist(link_vars)) {
                parent_lvl <- NULL
                for (lvl_name in names(data_list)) {
                    if (is.data.frame(data_list[[lvl_name]]) && link %in% names(data_list[[lvl_name]])) {
                        if (!any(duplicated(data_list[[lvl_name]][[link]]))) {
                            parent_lvl <- lvl_name
                            break
                        }
                    }
                }
                if (!is.null(parent_lvl)) {
                    parent_keys <- data_list[[parent_lvl]][[link]]
                    for (lvl_name in names(data_list)) {
                        if (is.data.frame(data_list[[lvl_name]]) && lvl_name != parent_lvl && link %in% names(data_list[[lvl_name]])) {
                            child_keys <- data_list[[lvl_name]][[link]]
                            idx_array <- match(child_keys, parent_keys) - 1 # 0-indexed for python
                            for (p_var in names(data_list[[parent_lvl]])) {
                                if (p_var != link) {
                                    flat_data[[paste0("idx_", p_var)]] <- idx_array
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    # Ensure all data passed to Python/JAX is numeric; drop character/factor columns (e.g. id_col)
    for (col in names(flat_data)) {
        if (!is.numeric(flat_data[[col]]) && !is.logical(flat_data[[col]])) {
            flat_data[[col]] <- NULL
        }
    }
    
    return(flat_data)
}
