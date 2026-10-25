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

  # Scale-aware cross-scale dispatch
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
    synth_iter <- max(n.iter, 1000)

    return(run_crossscale_dsep(
      i = i,
      test_eq = test_eq,
      cs_info = cs_info,
      original_data = original_data,
      hierarchical_info = hierarchical_info,
      structure = structure,
      family = family,
      engine = engine,
      n.iter = synth_iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      n.adapt = n.adapt,
      n.chains = n.chains,
      quiet = quiet
    ))
  }

  test_data <- original_data

  if (!is.null(hierarchical_info)) {
    test_vars <- all.vars(test_eq)

    # Track categorical parent variables whose dummy columns must be recreated after merging levels
    dummies_to_create <- list()
    cat_vars <- attr(original_data, "categorical_vars")

    if (!is.null(cat_vars)) {
      for (parent_var in names(cat_vars)) {
        if (parent_var %in% test_vars) {
          dummies_to_create[[parent_var]] <- cat_vars[[parent_var]]$dummies
        }
      }
    }

    test_data <- get_data_for_variables(
      test_vars,
      original_data,
      hierarchical_info$levels,
      hierarchical_info$hierarchy,
      hierarchical_info$link_vars,
      equations = equations,
      latent = latent
    )

    # Recreate dummy variables in merged test_data if missing
    for (parent_var in names(dummies_to_create)) {
      if (parent_var %in% names(test_data)) {
        dummies <- dummies_to_create[[parent_var]]
        levels_map <- cat_vars[[parent_var]]$levels
        if (is.null(levels_map)) next

        vals <- test_data[[parent_var]]
        type <- cat_vars[[parent_var]]$type

        if (!is.null(type) && type == "ordered") {
          c_mat <- cat_vars[[parent_var]]$contrasts
          match_idx <- if (is.numeric(vals)) vals else match(vals, levels_map)

          if (any(is.na(match_idx))) {
            warning(sprintf(
              "Categorical level mismatch: %d value(s) in '%s' not found in expected levels.",
              sum(is.na(match_idx)), deparse(substitute(vals))
            ))
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
            if (is.null(test_data[[expected_dummy]]) || all(is.na(test_data[[expected_dummy]]))) {
              is_match <- if (is.numeric(vals)) vals == k else vals == levels_map[k]
              test_data[[expected_dummy]] <- as.integer(is_match)
            }
          }
        }
      }
    }

    if (!is.null(cat_vars)) {
      attr(test_data, "categorical_vars") <- cat_vars
    }

    if (!quiet) {
      message("  Hierarchical data: using ", nrow(test_data), " observations for this test")
    }
  }

  # Expand d-separation equations via extension hook (e.g., detection sub-models)
  dsep_equations <- dsep_equations_hook(
    family,
    equations,
    list(test_eq),
    test_eq = test_eq
  )

  # Filter family to responses present in dsep_equations
  all_dsep_responses <- unique(unlist(lapply(dsep_equations, function(eq) as.character(eq)[2])))
  all_dsep_responses_clean <- unique(c(
    all_dsep_responses,
    sub("^p_", "", all_dsep_responses),
    sub("^psi_", "", all_dsep_responses),
    sub("^z_", "", all_dsep_responses)
  ))

  sub_family <- if (!is.null(family)) {
    family[names(family) %in% all_dsep_responses_clean]
  } else {
    NULL
  }

  # Auto-detect binomial family for binary (0/1) response variables
  test_resp <- as.character(test_eq)[2]
  if (!is.null(test_data[[test_resp]])) {
    vals <- if (is.list(test_data) && !is.data.frame(test_data)) {
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
        if (is.null(sub_family)) sub_family <- list()
        if (is.na(sub_family[test_resp])) sub_family[[test_resp]] <- "binomial"
      }
    }
  }

  if (length(sub_family) == 0) sub_family <- NULL

  dsep_data_to_pass <- if (!is.null(original_data)) original_data else test_data

  # Subsample observations if dataset exceeds dsep_max_obs
  if (!is.null(hierarchical_info)) {
    h_paths <- strsplit(multiscale, "\\s*;\\s*")[[1]]
    finest_lvls <- unique(sapply(h_paths, function(path) {
      lvls <- trimws(strsplit(path, "\\s*>\\s*")[[1]])
      lvls[length(lvls)]
    }))
    for (fl in finest_lvls) {
      if (fl %in% names(dsep_data_to_pass) && is.data.frame(dsep_data_to_pass[[fl]])) {
        n_obs <- nrow(dsep_data_to_pass[[fl]])
        if (n_obs > dsep_max_obs) {
          set.seed(i + 42)
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

  sub_structure <- dsep_tree_hook(structure, test_eq, hierarchical_info, levels)

  # Truncate hierarchy and aggregate finer-level variables when the test is locked to a coarse scale
  sub_multiscale <- multiscale
  sub_levels <- levels
  sub_link_vars <- link_vars

  if (!is.null(hierarchical_info)) {
    test_scale <- attr(test_eq, "scale")
    if (!is.null(test_scale)) {
      if (!is.null(dsep_data_to_pass) && is.list(dsep_data_to_pass) && !is.data.frame(dsep_data_to_pass)) {
        dsep_data_to_pass <- aggregate_multiscale_data(
          dsep_data_to_pass,
          levels,
          multiscale,
          test_scale
        )
        if (!quiet) message(sprintf("  (Resolution Locked: Aggregated test to '%s' scale)", test_scale))
      }

      # Truncate hierarchy paths at test_scale
      h_paths <- strsplit(multiscale, "\\s*;\\s*")[[1]]
      new_h_paths <- sapply(h_paths, function(path) {
        lvls <- trimws(strsplit(path, "\\s*>\\s*")[[1]])
        idx <- match(test_scale, lvls)
        if (!is.na(idx)) paste(lvls[1:idx], collapse = " > ") else NA
      })
      new_h_paths <- new_h_paths[!is.na(new_h_paths)]
      sub_multiscale <- paste(unique(new_h_paths), collapse = " ; ")

      all_keep_lvls <- unique(unlist(lapply(strsplit(sub_multiscale, "\\s*[>;]\\s*"), trimws)))
      sub_levels <- levels[names(levels) %in% all_keep_lvls]

      # Reassign variables aggregated from finer levels into test_scale
      t_depth <- get_level_depth(test_scale, multiscale)
      for (v in all.vars(test_eq)) {
        orig_lvl <- NULL
        for (lvl in names(levels)) {
          if (v %in% levels[[lvl]]) {
            orig_lvl <- lvl
            break
          }
        }

        if (!is.null(orig_lvl)) {
          v_depth <- get_level_depth(orig_lvl, multiscale)
          if (!is.na(v_depth) && !is.na(t_depth) && v_depth > t_depth) {
            sub_levels <- lapply(sub_levels, function(l) setdiff(l, v))
            if (test_scale %in% names(sub_levels)) {
              sub_levels[[test_scale]] <- unique(c(sub_levels[[test_scale]], v))
            }
          }
        }
      }

      if (!is.null(link_vars)) {
        sub_link_vars <- link_vars[names(link_vars) %in% names(sub_levels)]
      }
    }
  }

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
    variability = variability,
    engine = engine,
    nimble_samplers = nimble_samplers
  )

  fit <- do.call(because, call_args)

  param_map <- fit$parameter_map
  param_map$equation_index <- i

  list(
    samples = fit$samples,
    param_map = param_map,
    model = fit$model,
    test_index = i
  )
}
