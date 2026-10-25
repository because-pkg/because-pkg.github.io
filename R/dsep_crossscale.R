# Cross-Scale D-Separation via Coarse-Level Aggregation
#
# When a focal predictor lives at a *coarser* hierarchical level than the
# response (e.g., site-level or species-level predictor tested against
# observation-level response), this module aggregates the response and
# conditioning variables up to the predictor's hierarchical level and fits a
# native Bayesian model via because() (with any covariance structure applicable
# at that level).
#
# Primary functions:
#   - detect_crossscale_dsep()  detect coarser-level focal predictors
#   - run_crossscale_dsep()     aggregate + Bayesian fit at predictor's scale
#
# @keywords internal

#' Detect whether a d-sep test crosses hierarchical scales in the same chain
#'
#' A test is "cross-scale" when the focal predictor is at a coarser level than
#' the response within the SAME hierarchical chain (not orthogonal branches).
#' E.g., a site- or species-level variable predicting an obs-level response in
#' "site > survey > obs; species > obs".
#'
#' @param test_eq A formula carrying the attribute \code{test_var}.
#' @param hierarchical_info List with \code{$levels}, \code{$hierarchy},
#'   \code{$link_vars}.
#' @return A list. If cross-scale: \code{$is_crossscale = TRUE} plus
#'   \code{$response}, \code{$response_level}, \code{$predictor_level},
#'   \code{$test_var}. Otherwise \code{$is_crossscale = FALSE}.
#' @keywords internal
detect_crossscale_dsep <- function(test_eq, hierarchical_info) {
  if (is.null(hierarchical_info) ||
      is.null(hierarchical_info$levels) ||
      is.null(hierarchical_info$hierarchy)) {
    return(list(is_crossscale = FALSE))
  }

  test_var <- attr(test_eq, "test_var")
  if (is.null(test_var)) {
    rhs_vars <- all.vars(test_eq)[-1]
    if (length(rhs_vars) > 0L) test_var <- rhs_vars[1L]
  }
  if (is.null(test_var)) return(list(is_crossscale = FALSE))

  resp     <- as.character(test_eq)[2L]
  resp     <- sub("^psi_", "", resp)

  resp_lvl <- get_var_level_dsep(resp,     hierarchical_info)
  pred_lvl <- get_var_level_dsep(test_var, hierarchical_info)

  if (is.null(resp_lvl) || is.null(pred_lvl)) return(list(is_crossscale = FALSE))
  if (resp_lvl == pred_lvl)                   return(list(is_crossscale = FALSE))

  # Parse hierarchy chains (semicolon-separated, ">" = coarser > finer)
  paths <- strsplit(hierarchical_info$hierarchy, "\\s*;\\s*")[[1L]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1L]]))

  for (path in paths) {
    resp_idx <- match(resp_lvl, path)
    pred_idx <- match(pred_lvl, path)
    if (!is.na(resp_idx) && !is.na(pred_idx) && pred_idx < resp_idx) {
      # Predictor is coarser (smaller index) than response in the same chain
      return(list(
        is_crossscale   = TRUE,
        response        = resp,
        response_level  = resp_lvl,
        predictor_level = pred_lvl,
        test_var        = test_var
      ))
    }
  }
  list(is_crossscale = FALSE)
}

.cs_flatten_all <- function(original_data) {
  if (is.data.frame(original_data)) return(original_data)
  if (!is.list(original_data)) return(original_data)
  
  row_counts <- sapply(original_data, function(x) if (is.data.frame(x)) nrow(x) else 0)
  base_idx <- which.max(row_counts)
  if (length(base_idx) == 0) return(NULL)
  
  base_tbl <- original_data[[base_idx]]
  merged_tables <- c(base_idx)
  
  repeat {
    merged_in_this_pass <- FALSE
    for (i in seq_along(original_data)) {
      if (i %in% merged_tables) next
      tbl <- original_data[[i]]
      if (!is.data.frame(tbl)) next
      
      shared <- intersect(names(base_tbl), names(tbl))
      if (length(shared) > 0) {
        new_cols <- setdiff(names(tbl), names(base_tbl))
        if (length(new_cols) > 0) {
          base_tbl <- merge(base_tbl, tbl[, c(shared, new_cols), drop = FALSE], by = shared, all.x = TRUE)
        }
        merged_tables <- c(merged_tables, i)
        merged_in_this_pass <- TRUE
      }
    }
    # Count dataframes
    n_dfs <- sum(sapply(original_data, is.data.frame))
    if (!merged_in_this_pass || length(merged_tables) == n_dfs) break
  }
  
  base_tbl
}

#' Run a cross-scale d-sep test via coarse-level aggregation
#'
#' Implements the scale-aware d-sep testing strategy for hierarchical models:
#' \enumerate{
#'   \item Aggregate the response variable to the predictor's hierarchical level
#'         (link-scale mean for non-Gaussian families; arithmetic mean for Gaussian).
#'   \item Aggregate conditioning variables to the predictor's level.
#'   \item Filter \code{structure} to any covariance structures matching the
#'         predictor's level (via \code{get_order_labels_hook}).
#'   \item Fit the aggregated single-equation model natively via \code{because()}
#'         and return the posterior samples for the focal coefficient.
#' }
#'
#' @param i            Integer index of the d-sep test.
#' @param test_eq      Formula carrying \code{attr(., "test_var")}.
#' @param cs_info      List returned by \code{detect_crossscale_dsep()}.
#' @param original_data Hierarchical data list or flat data.frame.
#' @param hierarchical_info List with \code{$levels}, \code{$hierarchy},
#'   \code{$link_vars}.
#' @param structure    Named list of structure objects (e.g. spatial matrices,
#'   phylogenetic trees), as passed to \code{because()}.
#' @param family       Named character vector mapping response names to families.
#' @param engine       Bayesian engine (\code{"jags"}, \code{"nimble"}, or \code{"numpyro"}).
#' @param n.iter       MCMC iterations per chain (default 1000).
#' @param n.burnin     Burn-in iterations per chain (default 500).
#' @param n.thin       Thinning interval (default 1).
#' @param n.adapt      Adaptation iterations (default 500).
#' @param n.chains     Number of MCMC chains (default 3).
#' @param quiet        Suppress informational messages.
#' @return List with \code{$samples} (coda::mcmc.list), \code{$param_map},
#'   \code{$model}, \code{$test_index}.
#' @keywords internal
run_crossscale_dsep <- function(
  i,
  test_eq,
  cs_info,
  original_data,
  hierarchical_info,
  structure = NULL,
  family    = NULL,
  engine    = "numpyro",
  n.iter    = 1000L,
  n.burnin  = 500L,
  n.thin    = 1L,
  n.adapt   = 500L,
  n.chains  = 3L,
  quiet     = FALSE
) {
  resp       <- cs_info$response
  test_var   <- cs_info$test_var
  pred_level <- cs_info$predictor_level
  link_vars  <- hierarchical_info$link_vars
  group_col  <- link_vars[[pred_level]]

  if (is.null(group_col)) {
    stop(sprintf(
      "Cross-scale d-sep: no link_var defined for level '%s'. Cannot aggregate.",
      pred_level
    ))
  }

  # ---- 1. Use the globally flattened dataset ----
  flat <- .cs_flatten_all(original_data)

  # ---- 2. Aggregate response to the predictor's group level ----
  y_raw   <- flat[[resp]]
  grp     <- as.character(flat[[group_col]])
  fam_str <- if (!is.null(family) && resp %in% names(family)) family[[resp]] else "gaussian"

  agg_y <- if (fam_str %in% c("poisson", "negbin", "negbinomial", "zip", "zinb")) {
    tapply(log(pmax(y_raw, 0) + 0.5), grp, mean, na.rm = TRUE)
  } else if (fam_str %in% c("lognormal", "gamma", "exponential")) {
    tapply(log(pmax(y_raw, 1e-6)), grp, mean, na.rm = TRUE)
  } else if (fam_str %in% c("binomial", "bernoulli", "beta")) {
    tapply(y_raw, grp, function(x) {
      p <- mean(x, na.rm = TRUE)
      p <- pmax(0.01, pmin(0.99, p)) # Bound to prevent log(0)
      log(p / (1 - p)) # Logit transform
    })
  } else {
    tapply(y_raw, grp, mean, na.rm = TRUE)
  }
  group_ids <- names(agg_y)

  agg_df <- data.frame(
    stringsAsFactors = FALSE,
    setNames(list(group_ids, as.numeric(agg_y)), c(group_col, ".response"))
  )

  # ---- 3. Aggregate all conditioning variables to the group level ----
  # Exclude random effect grouping variables (which are just the link_vars)
  rhs_vars <- setdiff(all.vars(test_eq)[-1L], c(group_col, unlist(link_vars)))
  
  for (v in rhs_vars) {
    if (v %in% names(flat)) {
      v_raw <- flat[[v]]
      agg_v <- tapply(v_raw, grp, function(x) {
        if (is.numeric(x)) {
          mean(x, na.rm = TRUE)
        } else {
          # For categorical predictors, take the most frequent value (mode)
          tbl <- sort(table(x), decreasing = TRUE)
          if (length(tbl) > 0) names(tbl)[1] else NA
        }
      })
      # Ensure order matches agg_df
      agg_df[[v]] <- as.vector(agg_v[agg_df[[group_col]]])
    }
  }

  agg_df <- agg_df[stats::complete.cases(agg_df), , drop = FALSE]

  # ---- 4. Reformat aggregated dataset for flat Bayesian run ----
  avail_preds <- intersect(rhs_vars, names(agg_df))
  if (length(avail_preds) == 0L) {
    stop("No predictors available in aggregated data for cross-scale d-sep test.")
  }
  
  # Ensure the response has its original name for the because() call
  names(agg_df)[names(agg_df) == ".response"] <- resp
  
  # For structure label alignment in because(), rename the identifier column to '.raw_<group_col>'
  raw_group_col <- paste0(".raw_", group_col)
  names(agg_df)[names(agg_df) == group_col] <- raw_group_col
  
  # Keep only the necessary columns to avoid unrelated data issues
  keep_cols <- c(raw_group_col, resp, avail_preds)
  agg_df <- agg_df[, keep_cols, drop = FALSE]
  
  method_used <- sprintf("Bayesian aggregated model (%s)", if (tolower(engine) == "jags") "JAGS" else if (tolower(engine) == "nimble") "NIMBLE" else "NumPyro")
  if (!quiet) {
    message(sprintf(
      "  [Cross-scale d-sep] '%s' aggregated to %s level (n=%d) -> %s ... ",
      resp, pred_level, nrow(agg_df), method_used
    ), appendLF = FALSE)
    flush.console()
  }

  # Build the flat equation
  agg_formula <- stats::as.formula(
    paste(resp, "~", paste(avail_preds, collapse = " + "))
  )

  # ---- 5. Filter structure to only what applies to this level ----
  pruned_structure <- list()
  if (!is.null(structure) && is.list(structure)) {
    groups_present <- as.character(agg_df[[raw_group_col]])
    for (s_name in names(structure)) {
      s_obj <- structure[[s_name]]
      s_labels <- get_order_labels_hook(s_obj)
      if (is.null(s_labels)) {
        tree_obj <- get_tree_hook(s_obj)
        if (!is.null(tree_obj)) s_labels <- get_order_labels_hook(tree_obj)
      }
      if (!is.null(s_labels) && length(intersect(s_labels, groups_present)) > 0) {
        pruned_structure[[s_name]] <- s_obj
      }
    }
  }
  if (length(pruned_structure) == 0) pruned_structure <- NULL

  # ---- 6. Run native fully-Bayesian model via because() ----
  fit_bayesian <- tryCatch({
    because(
      equations = list(agg_formula),
      data = agg_df,
      structure = pruned_structure,  # Pass only applicable structures
      family = {
        if (!is.null(family)) {
          agg_fam <- as.list(family)
          agg_fam[[resp]] <- "gaussian"
          agg_fam
        } else {
          NULL
        }
      },
      dsep = FALSE,           # Just fit this single equation, don't generate more tests
      engine = engine,        # Use user-specified engine
      n.chains = n.chains,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      n.adapt = n.adapt,
      quiet = TRUE            # Suppress inner compilation messages
    )
  }, error = function(e) {
    warning(sprintf(
      "Cross-scale Bayesian d-sep test failed (%s); returning empty samples for test %d.", conditionMessage(e), i
    ))
    NULL
  })

  param_name <- paste0("beta_", resp, "_", test_var)
  
  if (!quiet) {
    if (!is.null(fit_bayesian)) {
      message("Done.")
    } else {
      message("Failed.")
    }
  }

  # ---- 7. Extract MCMC chains and map parameters ----
  if (!is.null(fit_bayesian) && !is.null(fit_bayesian$samples)) {
    if (inherits(fit_bayesian$samples, "mcmc.list")) {
      extracted_samples <- coda::mcmc.list(lapply(fit_bayesian$samples, function(x) {
        coda::mcmc(as.matrix(x)[, param_name, drop = FALSE])
      }))
    } else {
      extracted_samples <- fit_bayesian$samples[, param_name, drop = FALSE]
    }
  } else {
    extracted_samples <- lapply(seq_len(n.chains), function(ch) {
      mat <- matrix(NA_real_, nrow = n.iter, ncol = 1L, dimnames = list(NULL, param_name))
      coda::mcmc(mat)
    })
    extracted_samples <- coda::mcmc.list(extracted_samples)
  }

  # ---- 8. param_map (minimal columns consumed by downstream code) ----
  param_map <- data.frame(
    response       = resp,
    predictor      = test_var,
    parameter      = param_name,
    equation_index = i,
    type           = "beta",
    stringsAsFactors = FALSE
  )

  model_str <- sprintf(
    paste0(
      "# Cross-scale Bayesian aggregated d-sep test (test %d)\n",
      "# Response '%s' aggregated to '%s' level (n = %d)\n",
      "# Formula: %s ~ %s\n",
      "# Method: %s"
    ),
    i, resp, pred_level, nrow(agg_df),
    resp, paste(avail_preds, collapse = " + "),
    method_used
  )

  list(
    samples    = extracted_samples,
    param_map  = param_map,
    model      = model_str,
    test_index = i
  )
}
