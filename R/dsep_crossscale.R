
# ==========================================================================
# Cross-Scale D-Sep: PGLS routing
# ==========================================================================
# When the focal predictor lives at a *coarser* hierarchical level than the
# response (e.g., species-level Body_Mass_s tested against obs-level Abundance),
# a standard observation-level GLMM is non-identifiable: the group-level fixed
# effect and the group-level random effect compete for the same variance.
#
# To avoid MCMC convergence issues, we run the test at the predictor's own scale:
#   1. Aggregate the response to the predictor's level.
#   2. Include only conditioning variables at or above that level.
#   3. Fit GLS with the appropriate correlation structure (phylo / spatial).
# ==========================================================================

#' Detect whether a d-sep test crosses hierarchical scales in the same chain
#'
#' A test is "cross-scale" when the focal predictor is at a coarser level than
#' the response within the SAME hierarchical chain (not orthogonal branches).
#' E.g., species-level Body_Mass_s predicting obs-level Abundance in the chain
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

#' Run a cross-scale d-sep test via aggregation and PGLS/GLS
#'
#' Implements the scale-aware d-sep testing strategy for hierarchical models:
#' \enumerate{
#'   \item Aggregate the response variable to the predictor's hierarchical level
#'         (log-mean for Poisson; arithmetic mean for Gaussian).
#'   \item Select conditioning variables at or above the predictor's level.
#'   \item Fit a GLS with \code{ape::corPagel} (Pagel's lambda estimated by ML)
#'         if a phylogenetic tree matching the predictor's level is found in
#'         \code{structure}; otherwise fall back to ordinary OLS.
#'   \item Return synthetic MCMC samples (normal approximation to the GLS
#'         posterior of the focal coefficient) so that the result slots directly
#'         into because()'s existing downstream summary / plot code.
#' }
#'
#' @param i            Integer index of the d-sep test.
#' @param test_eq      Formula carrying \code{attr(., "test_var")}.
#' @param cs_info      List returned by \code{detect_crossscale_dsep()}.
#' @param original_data Hierarchical data list or flat data.frame.
#' @param hierarchical_info List with \code{$levels}, \code{$hierarchy},
#'   \code{$link_vars}.
#' @param structure    Named list of structure objects (phylo trees, spatial
#'   matrices), as passed to \code{because()}.
#' @param family       Named character vector mapping response names to families.
#' @param n.iter       Synthetic draws per MCMC chain (default 1000).
#' @param n.chains     Number of synthetic chains (default 3).
#' @param quiet        Suppress informational messages.
#' @return List with \code{$samples} (coda::mcmc.list), \code{$param_map},
#'   \code{$model}, \code{$test_index}.
#' @keywords internal
run_crossscale_dsep_pgls <- function(
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
    stop("No predictors available in aggregated data for cross-scale PGLS.")
  }
  
  # Ensure the response has its original name for the because() call
  names(agg_df)[names(agg_df) == ".response"] <- resp
  
  # For the phylogeny mapping in because(), if the data is flat, the row names or 
  # an identifier column matching the tree tips must be present.
  # The identifier column is group_col, so we rename it to the expected raw format.
  # because() expects '.raw_Species' for species mapping. We use group_col dynamically.
  raw_group_col <- paste0(".raw_", group_col)
  names(agg_df)[names(agg_df) == group_col] <- raw_group_col
  
  # Ensure we only keep the necessary columns to avoid unrelated data issues
  keep_cols <- c(raw_group_col, resp, avail_preds)
  agg_df <- agg_df[, keep_cols, drop = FALSE]
  
  method_used <- sprintf("Bayesian PGLS/GLS (%s)", if (tolower(engine) == "jags") "JAGS" else if (tolower(engine) == "nimble") "NIMBLE" else "NumPyro")
  if (!quiet) {
    message(sprintf(
      "  [Cross-scale d-sep] '%s' aggregated to %s level (n=%d) -> %s ... ",
      resp, pred_level, nrow(agg_df), method_used
    ), appendLF = FALSE)
    flush.console()
  }

  # Build the flat equation
  gls_formula <- stats::as.formula(
    paste(resp, "~", paste(avail_preds, collapse = " + "))
  )

  # ---- 5. Prune structure to only what applies to this level ----
  pruned_structure <- list()
  if (!is.null(structure) && is.list(structure)) {
    groups_present <- as.character(agg_df[[raw_group_col]])
    for (s_name in names(structure)) {
      s_obj <- structure[[s_name]]
      # For phylo trees
      if (inherits(s_obj, "phylo")) {
        common <- intersect(s_obj$tip.label, groups_present)
        if (length(common) > 0) {
          pruned_structure[[s_name]] <- s_obj
        }
      } 
      # For spatial distance matrices
      else if (is.matrix(s_obj)) {
        common <- intersect(rownames(s_obj), groups_present)
        if (length(common) > 0) {
          pruned_structure[[s_name]] <- s_obj
        }
      }
    }
  }
  if (length(pruned_structure) == 0) pruned_structure <- NULL

  # ---- 6. Run native fully-Bayesian model via because() ----
  fit_bayesian <- tryCatch({
    because(
      equations = list(gls_formula),
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
      "Cross-scale Bayesian PGLS failed (%s); returning empty samples for test %d.", conditionMessage(e), i
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

  # ---- 6. Extract MCMC chains and map parameters ----
  
  if (!is.null(fit_bayesian) && !is.null(fit_bayesian$samples)) {
    # Extract the chains for the focal predictor from the Bayesian run
    # because() names coefficients as `beta_Response_Predictor`
    if (inherits(fit_bayesian$samples, "mcmc.list")) {
      extracted_samples <- coda::mcmc.list(lapply(fit_bayesian$samples, function(x) {
        coda::mcmc(as.matrix(x)[, param_name, drop = FALSE])
      }))
    } else {
      extracted_samples <- fit_bayesian$samples[, param_name, drop = FALSE]
    }
  } else {
    # Fallback if the Bayesian model failed to compile/run
    extracted_samples <- lapply(seq_len(n.chains), function(ch) {
      mat <- matrix(NA_real_, nrow = n.iter, ncol = 1L, dimnames = list(NULL, param_name))
      coda::mcmc(mat)
    })
    extracted_samples <- coda::mcmc.list(extracted_samples)
  }

  # ---- 7. param_map (minimal columns consumed by downstream code) ----
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
      "# Cross-scale Bayesian PGLS d-sep test (test %d)\n",
      "# Response '%s' aggregated to '%s' level (n = %d)\n",
      "# Formula: %s ~ %s\n",
      "# Method: %s"
    ),
    i, resp, pred_level, nrow(agg_df),
    resp, paste(avail_preds, collapse = " + "),
    method_used
  )

  list(
    samples       = extracted_samples,
    param_map     = param_map,
    model         = model_str,
    test_index    = i,
    exact_p_value = NA_real_ # P-values are not standard in pure Bayesian outputs; credible intervals will be used
  )
}


# ---- Internal helpers (not exported) ------------------------------------

# Flatten a hierarchical data list to a single data.frame that contains both
# the response variable and the group column needed for aggregation.
.cs_flatten <- function(original_data, resp, group_col) {
  if (is.data.frame(original_data)) return(original_data)
  if (!is.list(original_data))
    stop("original_data must be a data.frame or a named list of data.frames.")

  # Find the table that holds the response
  base_tbl <- NULL
  for (tbl in original_data) {
    if (is.data.frame(tbl) && resp %in% names(tbl)) { base_tbl <- tbl; break }
  }
  if (is.null(base_tbl))
    stop(sprintf("Response variable '%s' not found in any data table.", resp))

  # If the group column is already present we are done
  if (group_col %in% names(base_tbl)) return(base_tbl)

  # Otherwise try to merge in the table that has group_col
  for (tbl in original_data) {
    if (!is.data.frame(tbl) || identical(tbl, base_tbl)) next
    if (group_col %in% names(tbl)) {
      shared   <- intersect(names(base_tbl), names(tbl))
      new_cols <- setdiff(names(tbl), names(base_tbl))
      if (length(shared) > 0L && length(new_cols) > 0L) {
        base_tbl <- merge(
          base_tbl,
          tbl[, c(shared, new_cols), drop = FALSE],
          by = shared, all.x = TRUE
        )
      }
      if (group_col %in% names(base_tbl)) break
    }
  }

  if (!group_col %in% names(base_tbl))
    stop(sprintf("Group column '%s' not found after attempting to flatten data.", group_col))

  base_tbl
}

# Find the data.frame (level table) that contains both group_col and test_var.
.cs_find_table <- function(original_data, group_col, test_var) {
  if (is.data.frame(original_data)) return(original_data)
  for (tbl in original_data) {
    if (is.data.frame(tbl) &&
        group_col %in% names(tbl) &&
        test_var  %in% names(tbl)) {
      return(tbl)
    }
  }
  stop(sprintf(
    "Cannot find a table containing both '%s' and '%s' for cross-scale PGLS.",
    group_col, test_var
  ))
}

# Build a phylogenetic correlation structure for nlme::gls.
# Returns list($cor, $agg_df) where $agg_df is row-aligned to tree tip order.
.cs_build_cor <- function(structure, agg_df, group_col) {
  if (is.null(structure) || !requireNamespace("ape", quietly = TRUE)) {
    return(list(cor = NULL, agg_df = agg_df))
  }

  groups <- as.character(agg_df[[group_col]])

  for (sn in names(structure)) {
    s_obj <- structure[[sn]]
    if (!inherits(s_obj, "phylo")) next

    tree   <- s_obj
    common <- intersect(tree$tip.label, groups)

    if (length(common) < 3L) {
      return(list(cor = NULL, agg_df = agg_df))
    }

    # Prune tree to the groups present in the data
    tree_p <- ape::drop.tip(tree, setdiff(tree$tip.label, groups))

    # Make ultrametric (required by corPagel)
    if (!ape::is.ultrametric(tree_p, tol = 1e-6)) {
      tree_p <- ape::compute.brlen(tree_p, method = "Grafen")
    }

    # Align agg_df rows to tree tip order
    rownames(agg_df) <- agg_df[[group_col]]
    agg_df_aligned   <- agg_df[tree_p$tip.label, , drop = FALSE]

    cor_struct <- ape::corPagel(1, phy = tree_p, form = stats::as.formula(paste("~", group_col)), fixed = FALSE)
    return(list(cor = cor_struct, agg_df = agg_df_aligned))
  }

  # No matching phylo structure found
  list(cor = NULL, agg_df = agg_df)
}
