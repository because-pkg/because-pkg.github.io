run_single_dsep_test_v2 <- function(
  i,
  test_eq,
  monitor_params,
  engine = "jags",
  nimble_samplers = NULL,
  quiet = FALSE,
  original_data = NULL,
  hierarchical_info = NULL,
  random_terms = list(),
  equations = list(),
  family = NULL,
  structure = NULL,
  levels = NULL,
  hierarchy = NULL,
  multiscale = NULL,
  link_vars = NULL,
  fix_residual_variance = NULL,
  latent = NULL,
  latent_method = "correlations",
  n.chains = 3,
  n.iter = 12500,
  n.burnin = 2500,
  n.thin = 10,
  n.adapt = 2500,
  ic_recompile = FALSE,
  random = NULL,
  id_col = NULL,
  variability = NULL,
  dsep_max_obs = 10000,
  aggregate_crossscale = NULL
) {
  if (!quiet) {
    message(paste("D-sep test eq:", deparse(test_eq)))
    message(paste("Test var attribute:", attr(test_eq, "test_var")))
  }

  # --- AUTOMATED SCALE-AWARE DISPATCH ---
  cs_info <- detect_crossscale_dsep(test_eq, hierarchical_info)
  if (!quiet) message("is_crossscale:", cs_info$is_crossscale)

  do_aggregate <- FALSE
  if (cs_info$is_crossscale) {
    if (!quiet) {
      message(sprintf(
        "  -> Cross-scale test detected: focal predictor '%s' at %s level, response '%s' at %s level",
        cs_info$test_var, cs_info$predictor_level, cs_info$response, cs_info$response_level
      ))
    }
    
    if (is.character(aggregate_crossscale) && length(aggregate_crossscale) == 1 && aggregate_crossscale == "all") {
      do_aggregate <- TRUE
    } else if (is.numeric(aggregate_crossscale) && i %in% aggregate_crossscale) {
      do_aggregate <- TRUE
    }
  }

  if (do_aggregate) {
    synth_iter <- n.iter
    if (synth_iter < 1000) synth_iter <- 1000
    
    return(run_crossscale_dsep_pgls(
      i = i,
      test_eq = test_eq,
      cs_info = cs_info,
      original_data = original_data,
      hierarchical_info = hierarchical_info,
      structure = structure,
      family = family,
      engine = engine,  # Pass the engine choice down!
      n.iter = synth_iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      n.adapt = n.adapt,
      n.chains = n.chains,
      quiet = quiet
    ))
  }
  # ----------------------------------------

  # Select appropriate dataset for this test
  test_data <- original_data

  if (!is.null(hierarchical_info)) {
    # Extract variables from this test equation
    test_vars <- all.vars(test_eq)
    
    # [FIX 2026-04-13] we NO LONGER strip link_vars from test_vars.
    # infer_variable_level() now has a fallback to search data columns,
    # and get_data_for_variables() needs these IDs to correctly assemble/join levels.
    
    # [FIX] Add random effect grouping variables to test_vars
    # Otherwise get_data_for_variables removes them, causing "Unknown variable N_SiteID"
    # [FIX] categorical dummy variables
    # We need their PARENT variables to fetch the data, then recreate dummies manually.
    dummies_to_create <- list()
    cat_vars <- NULL

    if (!is.null(attr(original_data, "categorical_vars"))) {
      cat_vars <- attr(original_data, "categorical_vars")

      # Check all cat vars to see if their dummies are needed (or if parent is needed)
      current_vars <- test_vars

      for (parent_var in names(cat_vars)) {
        dummies <- cat_vars[[parent_var]]$dummies

        # If parent variable is in the test variables, we need to ensure we can recreate expected dummies
        if (parent_var %in% current_vars) {
          dummies_to_create[[parent_var]] <- dummies
        }
      }
    }

    # Get appropriate dataset for these variables (dummies NOT included in request)
    test_data <- get_data_for_variables(
      test_vars,
      original_data,
      hierarchical_info$levels,
      hierarchical_info$hierarchy,
      hierarchical_info$link_vars,
      equations = equations,
      latent = latent
    )
    
    # [FIX] Preserve categorical metadata for K detection in sub-models
    if (!is.null(cat_vars)) {
      attr(test_data, "categorical_vars") <- cat_vars
    }

    # [FIX] Recreate dummy variables in test_data
    for (parent_var in names(dummies_to_create)) {
      if (parent_var %in% names(test_data)) {
        dummies <- dummies_to_create[[parent_var]]
        vals <- test_data[[parent_var]]

        # Recreate dummies: sex_m = as.integer(sex == 2) etc.
        cat_vars_attr <- attr(original_data, "categorical_vars")
        levels_map <- if (
          !is.null(cat_vars_attr) && parent_var %in% names(cat_vars_attr)
        ) {
          cat_vars_attr[[parent_var]]$levels
        } else {
          NULL
        }
        if (is.null(levels_map)) {
          next
        }

        # Create all dummies for this parent
        type <- cat_vars[[parent_var]]$type
        
        if (!is.null(type) && type == "ordered") {
          c_mat <- cat_vars[[parent_var]]$contrasts
          vals <- test_data[[parent_var]]
          
          # Handle values which might be numeric indices or strings
          match_idx <- if (is.numeric(vals)) {
            vals
          } else {
            match(vals, levels_map)
          }
          
          if (any(is.na(match_idx))) {
            warning(sprintf("Categorical level mismatch: %d value(s) in '%s' not found in expected levels.",
                            sum(is.na(match_idx)), deparse(substitute(vals))))
          }

          for (k in seq_along(dummies)) {
            expected_dummy <- dummies[k]
            if (is.null(test_data[[expected_dummy]]) || all(is.na(test_data[[expected_dummy]]))) {
              test_data[[expected_dummy]] <- c_mat[match_idx, k]
            }
          }
        } else {
          for (k in 2:length(levels_map)) {
            expected_dummy <- dummies[k - 1]
  
            # Check if it was extracted properly. If not, recreate it!
            # We only recreate if it's missing or NA
            if (
              is.null(test_data[[expected_dummy]]) ||
                all(is.na(test_data[[expected_dummy]]))
            ) {
              vals <- test_data[[parent_var]]
  
              # Robust check matching preprocess_categorical_vars logic
              is_match <- if (is.numeric(vals)) {
                vals == k
              } else {
                vals == levels_map[k]
              }
  
              test_data[[expected_dummy]] <- as.integer(is_match)
            }
          }
        }
      }
    }

    # [FIX] Restore categorical_vars attribute dropped by merge/get_data
    if (!is.null(attr(original_data, "categorical_vars"))) {
      attr(test_data, "categorical_vars") <- attr(
        original_data,
        "categorical_vars"
      )
    }

    if (!quiet) {
      message(
        "  Hierarchical data: using ",
        nrow(test_data),
        " observations for this test"
      )
    }
  }

  # [REVISED] Populate dsep_equations FIRST, then filter sub_family/sub_variability
  dsep_equations <- list(test_eq)

  test_eq_vars <- all.vars(test_eq)
  # [NEW 2025-12-21] For occupancy models, we must include supporting equations
  # (detection models and models for latent predictors)
  # [Fixed duplicate line]

  # Extension Hook: Expand d-separation equations (e.g. detection models)
  dsep_equations <- dsep_equations_hook(
    family,
    equations,
    dsep_equations,
    test_eq = test_eq
  )

  # Now filter variability/family based on ALL variables in dsep_equations
  # [Manual Fix: variability needs to be handled if it's there, but here we only have family]
  all_dsep_vars <- unique(unlist(lapply(dsep_equations, all.vars)))
  all_dsep_responses <- unique(unlist(lapply(dsep_equations, function(eq) as.character(eq)[2])))
  
  all_dsep_vars_clean <- unique(c(
    all_dsep_vars,
    sub("^p_", "", all_dsep_vars),
    sub("^psi_", "", all_dsep_vars),
    sub("^z_", "", all_dsep_vars)
  ))
  
  all_dsep_responses_clean <- unique(c(
    all_dsep_responses,
    sub("^p_", "", all_dsep_responses),
    sub("^psi_", "", all_dsep_responses),
    sub("^z_", "", all_dsep_responses)
  ))

  # Note: sub_variability logic removed for brevity if not strictly needed in this context
  # but we'll try to extract what we can from arguments
  sub_family <- if (!is.null(family)) {
    family[names(family) %in% all_dsep_responses_clean]
  } else {
    NULL
  }

  # [User Request] Auto-detect binomial for binary response
  test_resp <- as.character(test_eq)[2]
  if (!is.null(test_data[[test_resp]])) {
    # Try to find target column in test_data (which might be a list or df)
    vals <- if (is.list(test_data) && !is.data.frame(test_data)) {
      # find which element contains it
      found_vals <- NULL
      for (lvl in names(test_data)) {
        if (test_resp %in% names(test_data[[lvl]])) {
          found_vals <- test_data[[lvl]][[test_resp]]
          break
        }
      }
      found_vals
    } else {
      test_data[[test_resp]]
    }

    if (!is.null(vals)) {
      u_vals <- unique(na.omit(vals))
      if (length(u_vals) <= 2 && all(u_vals %in% c(0, 1))) {
        if (is.null(sub_family)) {
          sub_family <- list()
        }
        if (is.na(sub_family[test_resp])) {
          sub_family[[test_resp]] <- "binomial"
        }
      }
    }
  }

  if (length(sub_family) == 0) {
    sub_family <- NULL
  }

  # Choose what data to pass to the dsep sub-fit:
  dsep_data_to_pass <- if (!is.null(original_data)) {
    original_data
  } else {
    test_data
  }

  # [OPTIMIZATION] Down-sample for d-sep if data is massive (Diagnostic Performance)
  if (!is.null(hierarchical_info)) {
    # Down-sample the finest grain(s) in the list
    h_paths <- strsplit(multiscale, "\\s*;\\s*")[[1]]
    finest_lvls <- unique(sapply(h_paths, function(path) {
      lvls <- trimws(strsplit(path, "\\s*>\\s*")[[1]])
      lvls[length(lvls)]
    }))
    for (fl in finest_lvls) {
      if (fl %in% names(dsep_data_to_pass) && is.data.frame(dsep_data_to_pass[[fl]])) {
        n_obs <- nrow(dsep_data_to_pass[[fl]])
        if (n_obs > dsep_max_obs) {
          set.seed(i + 42) # Reproducible diagnostic sample
          dsep_data_to_pass[[fl]] <- dsep_data_to_pass[[fl]][sample(1:n_obs, dsep_max_obs), ]
          if (!quiet) message(sprintf("  (Optimized d-sep: subsampled '%s' scale to %d rows)", fl, dsep_max_obs))
        }
      }
    }
  } else if (is.data.frame(dsep_data_to_pass)) {
    n_obs <- nrow(dsep_data_to_pass)
    if (n_obs > dsep_max_obs) {
      set.seed(i + 42)
      dsep_data_to_pass <- dsep_data_to_pass[sample(1:n_obs, dsep_max_obs), ]
      if (!quiet) message(sprintf("  (Optimized d-sep: subsampled to %d rows)", dsep_max_obs))
    }
  }

  # Extension Hook: Filter structure
  sub_structure <- dsep_tree_hook(structure, test_eq, hierarchical_info, levels)

  # [NEW] Multiscale Resolution Locking: Truncate hierarchy for d-sep sub-models.
  # If the test is performed at a coarse level (e.g. Year), the sub-model should
  # not expect or look for Individual-level data/loops.
  sub_multiscale <- multiscale
  sub_levels <- levels
  if (!is.null(hierarchical_info)) {
    test_scale <- attr(test_eq, "scale")
    if (!is.null(test_scale)) {
      # [RESOLUTION LOCKING] Aggregate finer variables up to the test scale.
      # This prevents DOF inflation by ensuring we don't treat fine-grained
      # observations as independent points when the test is at a coarser scale.
      if (!is.null(dsep_data_to_pass) && is.list(dsep_data_to_pass) && !is.data.frame(dsep_data_to_pass)) {
          dsep_data_to_pass <- aggregate_multiscale_data(
              dsep_data_to_pass,
              sub_levels,
              multiscale, # Use original multiscale for depth checking
              test_scale
          )
          if (!quiet) message(sprintf("  (Resolution Locked: Aggregated test to '%s' scale)", test_scale))
      }

      # Truncate hierarchy string to stop at test_scale
      h_paths <- strsplit(multiscale, "\\s*;\\s*")[[1]]
      new_h_paths <- sapply(h_paths, function(path) {
        lvls <- trimws(strsplit(path, "\\s*>\\s*")[[1]])
        idx <- match(test_scale, lvls)
        if (!is.na(idx)) {
          paste(lvls[1:idx], collapse = " > ")
        } else {
          NA # Drop the path completely if it doesn't contain the test scale!
        }
      })
      new_h_paths <- new_h_paths[!is.na(new_h_paths)]
      sub_multiscale <- paste(unique(new_h_paths), collapse = " ; ")

      # Filter levels to only include those present in the truncated hierarchy
      all_keep_lvls <- unique(unlist(lapply(strsplit(sub_multiscale, "\\s*[>;]\\s*"), trimws)))
      sub_levels <- levels[names(levels) %in% all_keep_lvls]

      # [AGGREGATION] Actually aggregate the data to the test scale.
      # This ensures that variables from finer levels are available in the
      # dataframe of the test scale, so validate_hierarchical_data doesn't fail.
      sub_data_list <- aggregate_multiscale_data(original_data, levels, multiscale, test_scale)
      dsep_data_to_pass <- sub_data_list
      
      # Now re-assign variables from finer levels into the coarser test scale level
      # ONLY if they were actually aggregated (i.e., their original depth was > test_scale depth).
      # Coarser predictors are broadcasted downwards via indices and should remain in their original levels.
      t_depth <- get_level_depth(test_scale, multiscale)
      vars_in_test <- all.vars(test_eq)
      for (v in vars_in_test) {
        # Find original level of this variable
        orig_lvl <- NULL
        for (lvl in names(levels)) {
          if (v %in% levels[[lvl]]) {
            orig_lvl <- lvl
            break
          }
        }
        
        if (!is.null(orig_lvl)) {
          v_depth <- get_level_depth(orig_lvl, multiscale)
          # If the variable was aggregated upwards from a finer level
          if (!is.na(v_depth) && !is.na(t_depth) && v_depth > t_depth) {
            # Remove from any existing levels
            sub_levels <- lapply(sub_levels, function(l) setdiff(l, v))
            # Add to the effective test scale level
            if (test_scale %in% names(sub_levels)) {
              sub_levels[[test_scale]] <- unique(c(sub_levels[[test_scale]], v))
            }
          }
        }
      }

      # [FIX] Subset link_vars to only include levels present in the truncated hierarchy
      if (!is.null(link_vars)) {
        sub_link_vars <- link_vars[names(link_vars) %in% names(sub_levels)]
      } else {
        sub_link_vars <- NULL
      }
    }
  }

  if (!exists("sub_link_vars", inherits = FALSE)) sub_link_vars <- link_vars

  # Call because recursively
  # Use do.call and filtering to handle potential version conflicts on worker nodes
  bec_args <- names(formals(because))

  call_args <- list(
    data = dsep_data_to_pass,
    structure = sub_structure,
    equations = dsep_equations,
    monitor = monitor_params,
    n.chains = n.chains,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = n.thin,
    DIC = FALSE,
    WAIC = FALSE,
    n.adapt = n.adapt,
    quiet = quiet,
    dsep = FALSE,
    family = sub_family,
    fix_residual_variance = fix_residual_variance,
    latent = latent,
    latent_method = latent_method,
    parallel = FALSE,
    n.cores = 1,
    cl = NULL,
    ic_recompile = ic_recompile,
    random = random,
    levels = sub_levels,
    hierarchy = sub_multiscale,
    multiscale = sub_multiscale,
    link_vars = sub_link_vars,
    hierarchical_info = hierarchical_info,
    structure_multi = hierarchical_info$structure_multi,
    structure_levels = hierarchical_info$structure_levels,
    id_col = id_col,
    variability = variability
  )

  # Only add NIMBLE arguments if the version of because() on this node supports them
  if ("engine" %in% bec_args) {
    call_args$engine <- engine
    call_args$nimble_samplers <- nimble_samplers
  }

  fit <- do.call(because, call_args)

  # Extract samples, map, and model
  samples <- fit$samples
  model_string <- fit$model
  
  # DEBUG: Print model string if it failed (it won't get here if it failed in do.call)
  # Actually, let's print it BEFORE do.call by looking at what because() would generate.
  # Or just let it fail and I'll add a print in because() itself if quiet=FALSE.

  param_map <- fit$parameter_map

  # Update equation index in parameter map to match the d-sep test index
  param_map$equation_index <- i

  list(
    samples = samples,
    param_map = param_map,
    model = model_string,
    test_index = i
  )
}
