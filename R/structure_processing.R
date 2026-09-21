#' Process covariance structures for because models
#'
#' Validates and formats user-supplied covariance structures (spatial, phylogenetic,
#' etc.) into the data and metadata needed by JAGS/NIMBLE/NumPyro.
#' @keywords internal
process_because_structures <- function(
  structure, structures, structure_obj, data, hierarchical_info, is_hierarchical,
  equations, family, family_obj, levels, multiscale, latent, quiet,
  engine = "jags", original_raw_data_before_preprocess = NULL, row_ids = NULL
) {
  # --- Structure Processing ---
  structures <- list()
  is_multiple <- FALSE
  N <- NULL

  # 1. Normalize Input to List
  if (is.null(structure)) {
    # Independent model - no structures to process
  } else if (is.matrix(structure)) {
    structures[["custom"]] <- structure
  } else if (is.list(structure) && !inherits(structure, "list")) {
    # It's an S3 object with list base - use class name
    class_name <- class(structure)[1]
    # Check if it's a multi-object type (contains multiple items)
    # Defense: Ignore list-based S3 objects that represent a single entity (like 'phylo' or 'spatial_knn')
    if (length(structure) > 1 && is.null(names(structure))) {
      is_multiple <- TRUE
      N_trees <- length(structure)
    }
    structures[[class_name]] <- structure
  } else if (
    is.list(structure) && (is.null(class(structure)) || identical(class(structure), "list"))
  ) {
    # Plain list of structures
    structures <- structure
    # Check for multi-objects in the list
    for (s in structures) {
      if (is.list(s) && length(s) > 1 && !is.matrix(s) && !inherits(s, "phylo") && !inherits(s, "because_structure")) {
        # Generic list with multiple items (likely replicates or multiPhylo)
        is_multiple <- TRUE
        if (!exists("N_trees")) N_trees <- length(structures)
      }
    }
  } else {
    # Any other S3 object (phylo, spatial_knn, etc.) - use class name
    class_name <- class(structure)[1]
    structures[[class_name]] <- structure
  }

  # Ensure N_trees is available if is_multiple is true
  if (is_multiple && !exists("N_trees")) {
     N_trees <- length(structures)
  }

  # Discover total N (number of observations) early
  if (is.null(N) && "N" %in% names(data)) {
    N <- if (is.list(data)) data$N[1] else data[["N"]][1]
  }
  if (is.null(N)) {
    # Try to get N from data - handle both data.frame and list cases
    if (is.data.frame(data) && nrow(data) > 0) {
      N <- nrow(data)
    } else if (is.list(data) && length(data) > 0) {
      # For hierarchical lists, N should be the row count of the FINEST grain level.
      first_obj <- data[[1]]
      if (is.data.frame(first_obj)) {
        N <- nrow(first_obj)
      } else if (is.vector(first_obj) || is.factor(first_obj)) {
        N <- length(first_obj)
      } else if (is.matrix(first_obj) || is.array(first_obj)) {
        N <- nrow(first_obj)
      }
    }
  }
  # [CLEANUP] Only include N if we are not in a complex hierarchical context where level-specific Ns take priority.
  # This silences the "Unused variable 'N' in data" warning in JAGS.
  if (is.null(hierarchical_info) && !"N" %in% names(data)) {
    data$N <- N
  }

  # [URGENT FIX] Pre-calculate level-specific counts so structure auto-detection works!
  # If we don't do this here, we can't match e.g. a 50x50 matrix to the 50-species level.
  if (!is.null(hierarchical_info)) {
    for (lvl_name in names(hierarchical_info$levels)) {
      z_name <- paste0("zeros_", lvl_name)
      if (is.null(data[[z_name]])) {
        data[[z_name]] <- rep(0, N)
      }
      n_name <- paste0("N_", lvl_name)
      if (is.null(data[[n_name]])) {
        if (is.data.frame(data) && !is.null(data[[lvl_name]])) {
          data[[n_name]] <- as.integer(length(unique(na.omit(data[[lvl_name]]))))[1]
        } else if (is.list(data) && !is.null(data[[lvl_name]])) {
          if (is.data.frame(data[[lvl_name]]) || is.matrix(data[[lvl_name]])) {
            data[[n_name]] <- as.integer(nrow(data[[lvl_name]]))[1]
          } else {
            data[[n_name]] <- as.integer(length(data[[lvl_name]]))[1]
          }
        } else {
          data[[n_name]] <- as.integer(N)[1]
        }
      }
      if (!is.null(data[[n_name]])) {
        data[[n_name]] <- as.integer(data[[n_name]])[1]
      }
      
      # Populate synonymous names (ID column names)
      lvl_vars <- hierarchical_info$levels[[lvl_name]]
      for (v_nm in lvl_vars) {
        v_title <- paste0(toupper(substring(v_nm, 1, 1)), substring(v_nm, 2))
        potential_names <- unique(c(v_nm, v_title, toupper(v_nm), paste0(v_nm, "ID"), paste0(v_title, "ID"), paste0(v_nm, "_id"), paste0(v_title, "_id")))
        for (p_nm in potential_names) {
            nn_name <- paste0("N_", p_nm)
            zz_name <- paste0("zeros_", p_nm)
            if (is.null(data[[nn_name]])) data[[nn_name]] <- data[[n_name]]
            if (is.null(data[[zz_name]])) data[[zz_name]] <- data[[z_name]]
        }
      }
    }
  }

  structure_names <- names(structures)
  if (is.null(structure_names) && length(structures) > 0) {
    structure_names <- paste0("Struct", seq_along(structures))
    names(structures) <- structure_names
  }

  # 2. Process Structures using S3 Generic
  if (length(structures) == 0) {
    # Independent Logic: Determine N from data
    if (is.null(N) || N == 0) {
      potential_objects <- Filter(
        function(x) is.vector(x) || is.factor(x) || is.matrix(x) || is.array(x),
        data
      )
      if (length(potential_objects) > 0) {
        obj <- potential_objects[[1]]
        N <- if (is.matrix(obj) || is.array(obj)) nrow(obj) else length(obj)
      }
    }
  } else {
    # --- HOTFIX: Inject raw string labels so prepare_structure_data can align tips ---
    # JAGS data is flattened to integers, which destroys tip labels. 
    # We temporarily inject them back into 'data' so they can be discovered.
    injected_raw_cols <- character(0)
    if (is_hierarchical && !is.null(hierarchical_info$structure_levels)) {
      for (s_name in names(hierarchical_info$structure_levels)) {
        lvl_name <- hierarchical_info$structure_levels[[s_name]]
        link_var <- if (!is.null(hierarchical_info$link_vars)) hierarchical_info$link_vars[[lvl_name]] else NULL
        if (!is.null(link_var) && is.data.frame(original_raw_data_before_preprocess[[lvl_name]])) {
           raw_col_name <- paste0(".raw_", link_var)
           data[[raw_col_name]] <- as.character(original_raw_data_before_preprocess[[lvl_name]][[link_var]])
           injected_raw_cols <- c(injected_raw_cols, raw_col_name)
        }
      }
    }

    structure_levels <- list()
    structure_multi <- list()
    # Use S3 Generic for Processing
    for (s_name in structure_names) {
      structure_obj <- structures[[s_name]]
      prep_res <- prepare_structure_data(structure_obj, data = data, optimize = TRUE, quiet = quiet, engine = engine, row_ids = row_ids)

      if (!is.null(prep_res$data_list)) {
        for (d_name in names(prep_res$data_list)) {
          custom_s_name <- get_structure_name_hook(structure_obj)
          prefixed_name <- if (d_name %in% c("Prec", "VCV", "multiVCV", custom_s_name)) paste0(d_name, "_", s_name) else d_name
          
          # --- HOTFIX: Prevent Python/JAGS crashes ---
          # prepare_structure_data might return character/factor vectors (e.g. aligned tip labels).
          # We MUST NOT let these overwrite the integer index arrays (e.g. data$Species)
          # nor be added to the data list, as NumPyro strictly requires numeric arrays.
          if (is.character(prep_res$data_list[[d_name]]) || is.factor(prep_res$data_list[[d_name]])) {
             next
          }
          
          data[[prefixed_name]] <- prep_res$data_list[[d_name]]
        }
      } else {
        structures[[s_name]] <- NULL
        structure_names <- setdiff(structure_names, s_name)
        next
      }

      current_N <- NULL
      is_this_one_multi <- FALSE
      for (obj in prep_res$data_list) {
        if (is.matrix(obj) && nrow(obj) == ncol(obj)) {
          current_N <- nrow(obj)
          break
        } else if (is.array(obj) && length(dim(obj)) == 3) {
          dims <- dim(obj)
          if (dims[1] == dims[2]) {
             current_N <- dims[1]
             is_this_one_multi <- TRUE
             if (is_multiple && !exists("N_trees")) N_trees <- dims[3]
             break
          } else if (dims[2] == dims[3]) {
             current_N <- dims[2]
             is_this_one_multi <- TRUE
             if (is_multiple && !exists("N_trees")) N_trees <- dims[1]
             break
          }
        }
      }

      if (is.null(current_N)) {
        n_attr <- attr(structure_obj, "n")
        if (is.numeric(n_attr)) current_N <- n_attr
      }

      # [FIX] Robust Level Matching
      # Prioritize name-based matching (e.g. if the structure in the tree list is named "Species")
      s_level <- NULL
      if (!is.null(hierarchical_info) && !is.null(current_N)) {
        # 1. Try Name Matching first
        if (s_name %in% names(hierarchical_info$levels)) {
            s_level <- s_name
        } else {
            # 2. Fallback to Dimension Matching
            for (lvl in names(hierarchical_info$levels)) {
              n_name <- paste0("N_", lvl)
              if (!is.null(data[[n_name]]) && data[[n_name]] == current_N) {
                s_level <- lvl
                break
              }
            }
        }
      }
      structure_levels[[s_name]] <- s_level
      structure_multi[[s_name]]  <- is_this_one_multi

      if (!is.null(current_N)) {
        if (is.null(N) || N == 0) {
          N <- current_N
        } else if (N != current_N && is.null(hierarchical_info)) {
          stop(paste("Dimension mismatch in structure:", s_name))
        }
      }

      if (!is.null(prep_res$structure_object)) {
        structures[[s_name]] <- prep_res$structure_object
      }
    }

    if (!is.null(hierarchical_info)) {
      # Merge: prefer the already name-matched auto-detect result for each s_name,
      # only overwriting with the newly computed (dimension-based) value when the
      # auto-detect didn't produce a result for that structure.
      existing_sl <- hierarchical_info$structure_levels
      for (s_name in names(structure_levels)) {
        # Only override if auto-detect has no entry or if dimension match is more specific
        if (is.null(existing_sl[[s_name]])) {
          existing_sl[[s_name]] <- structure_levels[[s_name]]
        }
        # If both give a result, keep the existing (name-based) auto-detect result
      }
      hierarchical_info$structure_levels <- existing_sl
      hierarchical_info$structure_multi  <- structure_multi
    }
  }


  return(list(data = data, structures = structures))
}
