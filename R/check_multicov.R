#' Check Multi-Structure Covariance Convergence and Sensitivity Diagnostics
#'
#' Evaluates structure sampling, within-structure MCMC convergence (`Rhat_MCMC`),
#' within-group structure sensitivity, and between-group divergence (`Rhat_group`),
#' alongside an exact Law of Total Variance decomposition partitioning posterior
#' uncertainty into MCMC estimation variance, within-group structure variation, and
#' between-group conflict for any `because` model fitted with a collection or ensemble
#' of covariance structures (e.g., phylogenetic trees/backbones, spatial distance kernels,
#' genomic relationship matrices, or network graphs).
#'
#' @param model A `because` model object fitted with multiple candidate covariance structures.
#' @param plot Logical. If `TRUE` (default), displays a publication-ready 2-panel diagnostic
#'   figure featuring convergence/sensitivity forest plot (left) and 100% stacked uncertainty
#'   variance decomposition (right).
#' @param ... Additional arguments passed to plotting functions or metadata overrides
#'   (e.g., `group = c(...)` to override group/backbone labels).
#'
#' @return An invisible object of class `check_multicov` containing per-chain
#'   structure statistics, structure frequencies, parameter-level split
#'   `Rhat_MCMC` vs. `Rhat_group` diagnostics, and the exact uncertainty
#'   `variance_decomposition` table.
#'
#' @export
#' @importFrom graphics plot lines legend par grid barplot mtext abline points segments axis
#' @importFrom stats var sd quantile
#'
#' @examples
#' \dontrun{
#' # Fit with multiple candidate covariance structures
#' fit <- because(list(y ~ x), data = my_data, structure = my_cov_list, engine = "numpyro")
#' check_multicov(fit)
#' }
check_multicov <- function(model, plot = TRUE, ...) {
  if (!inherits(model, "because") && !is.list(model$samples)) {
    stop("Argument 'model' must be a 'because' model object containing '$samples'.")
  }

  if (is.null(model$samples) || length(model$samples) == 0) {
    stop("Model does not contain samples. Ensure the model was run with n.iter > 0.")
  }

  all_cols <- colnames(model$samples[[1]])
  n_chains <- length(model$samples)
  n_iters  <- nrow(model$samples[[1]])

  # 1. Locate structure/tree index column (e.g. K_struct_*, K_tree_*)
  k_col <- grep("^K_(struct|tree)_", all_cols, value = TRUE)
  if (length(k_col) == 0) {
    k_col <- grep("K_(struct|tree)", all_cols, value = TRUE)
  }

  dots <- list(...)

  if (length(k_col) == 0) {
    # Check if this model is an ensemble where each chain was mapped to a candidate structure
    n_structs_cand <- model$n_trees_total %||% model$n_structs_total
    if (is.null(n_structs_cand) && !is.null(model$structure) && (is.list(model$structure) || inherits(model$structure, "multiPhylo"))) {
      n_structs_cand <- length(model$structure)
    }
    if (is.null(n_structs_cand) && !is.null(model$structures)) {
      coll_in_struct <- Filter(function(s) is.list(s) || inherits(s, "multiPhylo"), model$structures)
      if (length(coll_in_struct) > 0) n_structs_cand <- length(coll_in_struct[[1]])
    }

    if (!is.null(n_structs_cand) && n_structs_cand > 1L && (n_chains %% n_structs_cand == 0L)) {
      cpt <- as.integer(n_chains / n_structs_cand)
      assigned_structs <- rep(seq_len(n_structs_cand), each = cpt)
      chain_data <- lapply(seq_len(n_chains), function(ch) {
        rep(as.integer(assigned_structs[ch]), n_iters)
      })
      k_col <- "K_struct_ensemble"
    } else {
      stop("No structure index parameter (e.g. 'K_struct_*' or 'K_tree_*') found in model samples.\n",
           "Was this model fitted with multiple candidate covariance structures?")
    }
  } else {
    k_col <- k_col[1]
    chain_data <- lapply(seq_len(n_chains), function(ch) {
      as.integer(round(model$samples[[ch]][, k_col]))
    })
  }

  # Detect whether this model represents phylogenetic structures
  is_phylo <- isTRUE(dots$is_phylo) ||
    grepl("tree|phylo", k_col, ignore.case = TRUE) ||
    inherits(model$structure, "multiPhylo") ||
    inherits(model$structure, "phylo") ||
    (is.list(model$structure) && length(model$structure) > 0 &&
       (inherits(model$structure[[1]], "phylo") || inherits(model$structure[[1]], "multiPhylo") ||
        (is.list(model$structure[[1]]) && length(model$structure[[1]]) > 0 && inherits(model$structure[[1]][[1]], "phylo")))) ||
    (!is.null(model$structures) && any(vapply(model$structures, function(s) {
      inherits(s, "multiPhylo") || inherits(s, "phylo") ||
        (is.list(s) && length(s) > 0 && inherits(s[[1]], "phylo"))
    }, logical(1))))

  s_term   <- dots$s_term %||% (if (is_phylo) "Tree" else "Structure")
  s_term_l <- tolower(s_term)
  g_term   <- dots$g_term %||% dots$group_name %||% dots$source_name %||% (if (is_phylo) "Backbone" else "Group")
  g_term_l <- tolower(g_term)

  # 2. Compute per-chain statistics
  stats_list <- lapply(seq_len(n_chains), function(ch) {
    structs <- chain_data[[ch]]
    unique_s <- sort(unique(structs))
    n_jumps <- sum(diff(structs) != 0)
    jump_rate <- n_jumps / max(1, length(structs) - 1)

    data.frame(
      Chain = ch,
      Unique_Structures = length(unique_s),
      Visited_Structures = paste(unique_s, collapse = ", "),
      Jumps = n_jumps,
      Jump_Rate_Pct = round(jump_rate * 100, 2),
      stringsAsFactors = FALSE
    )
  })
  chain_stats <- do.call(rbind, stats_list)

  # 3. Across-chain summary
  all_visited <- sort(unique(unlist(chain_data)))
  total_unique <- length(all_visited)
  total_jumps  <- sum(chain_stats$Jumps)

  # Cross-tabulation: Structure frequencies per chain
  crosstab_mat <- matrix(0L, nrow = length(all_visited), ncol = n_chains,
                         dimnames = list(paste0(s_term, "_", all_visited),
                                         paste0("Chain_", seq_len(n_chains))))
  for (ch in seq_len(n_chains)) {
    tbl <- table(factor(chain_data[[ch]], levels = all_visited))
    crosstab_mat[, ch] <- as.integer(tbl)
  }

  struct_totals <- rowSums(crosstab_mat)
  struct_weights_pct <- round(100 * struct_totals / sum(struct_totals), 2)

  # 4. Determine Mode: Deterministic Ensemble vs. Categorical Gibbs Jumping
  call_str <- if (!is.null(model$call)) paste(deparse(model$call), collapse = " ") else ""
  explicit_dcat <- identical(model$tree_method, "dcat") || identical(model$structure_method, "dcat") || grepl("dcat", call_str)
  explicit_ensemble <- identical(model$tree_method, "ensemble") || identical(model$structure_method, "ensemble")
  is_balanced_ensemble <- !explicit_dcat && (
    explicit_ensemble ||
      (total_jumps == 0 && total_unique > 1 &&
         max(struct_weights_pct) - min(struct_weights_pct) < 1.0 &&
         (is.null(model$n_trees_total) || total_unique == model$n_trees_total))
  )

  # 5. Discover Group / Backbone Metadata
  group <- NULL
  dots <- list(...)
  if (!is.null(dots$group)) {
    group <- as.character(dots$group)
  } else if (!is.null(dots$backbone)) {
    group <- as.character(dots$backbone)
  } else if (!is.null(model$structure) && !is.null(attr(model$structure, "group"))) {
    group <- as.character(attr(model$structure, "group"))
  } else if (!is.null(model$structure) && !is.null(attr(model$structure, "backbone"))) {
    group <- as.character(attr(model$structure, "backbone"))
  } else if (!is.null(model$input$structure) && !is.null(attr(model$input$structure, "group"))) {
    group <- as.character(attr(model$input$structure, "group"))
  } else if (!is.null(model$input$structure) && !is.null(attr(model$input$structure, "backbone"))) {
    group <- as.character(attr(model$input$structure, "backbone"))
  } else if (!is.null(model$structures)) {
    for (s in model$structures) {
      if (!is.null(attr(s, "group"))) {
        group <- as.character(attr(s, "group"))
        break
      }
      if (!is.null(attr(s, "backbone"))) {
        group <- as.character(attr(s, "backbone"))
        break
      }
    }
  }

  unique_groups <- if (!is.null(group)) unique(group) else character(0)
  has_groups <- length(unique_groups) > 1

  chain_primary_struct <- vapply(chain_data, function(x) x[1], integer(1))
  structs_with_multi_chains <- all_visited[vapply(all_visited, function(sk) {
    sum(chain_primary_struct == sk & chain_stats$Unique_Structures == 1) >= 2
  }, logical(1))]

  chain_group_vec <- if (has_groups) {
    vapply(chain_primary_struct, function(sk) {
      if (sk <= length(group)) group[sk] else paste0("Group_", sk)
    }, character(1))
  } else {
    rep("Default", n_chains)
  }

  # 6. Statistical helper functions
  calc_rhat <- function(mat_list) {
    m <- length(mat_list)
    if (m < 2) return(NA_real_)
    n <- length(mat_list[[1]])
    if (n < 2) return(NA_real_)
    c_means <- vapply(mat_list, mean, numeric(1), na.rm = TRUE)
    c_vars  <- vapply(mat_list, stats::var, numeric(1), na.rm = TRUE)
    W <- mean(c_vars, na.rm = TRUE)
    B <- n * stats::var(c_means, na.rm = TRUE)
    if (is.na(W) || W <= 1e-12) return(1.0)
    var_plus <- ((n - 1) / n) * W + (1 / n) * B
    sqrt(var_plus / W)
  }

  calc_group_rhat <- function(group_draws_list) {
    m <- length(group_draws_list)
    if (m < 2) return(NA_real_)
    g_means <- vapply(group_draws_list, mean, numeric(1), na.rm = TRUE)
    g_vars  <- vapply(group_draws_list, stats::var, numeric(1), na.rm = TRUE)
    n_per_g <- vapply(group_draws_list, length, integer(1))
    n_bar   <- mean(n_per_g)
    if (n_bar < 2) return(NA_real_)
    W <- mean(g_vars, na.rm = TRUE)
    B <- n_bar * stats::var(g_means, na.rm = TRUE)
    if (is.na(W) || W <= 1e-12) return(1.0)
    var_plus <- ((n_bar - 1) / n_bar) * W + (1 / n_bar) * B
    sqrt(var_plus / W)
  }

  # 7. Parameter Diagnostics
  struct_cols <- all_cols[
    grepl("^(alpha|beta|sigma|lambda|rho|psi|r_|cutpoint)", all_cols) &
      !grepl("^beta_.*(_det_|_1_times_)", all_cols)
  ]

  param_diag_df <- NULL
  if (length(struct_cols) > 0 && is_balanced_ensemble) {
    diag_rows <- lapply(struct_cols, function(p) {
      ch_vecs <- lapply(seq_len(n_chains), function(ch) as.numeric(model$samples[[ch]][, p]))
      all_draws <- unlist(ch_vecs)

      # -------------------------------------------------------------
      # Law of Total Variance (Sum of Squares Decomposition)
      # -------------------------------------------------------------
      grand_mean <- mean(all_draws, na.rm = TRUE)
      ss_total <- sum((all_draws - grand_mean)^2, na.rm = TRUE)

      draw_struct <- rep(chain_primary_struct, each = n_iters)
      struct_means <- tapply(all_draws, draw_struct, mean, na.rm = TRUE)
      draw_struct_mean <- struct_means[as.character(draw_struct)]
      ss_mcmc <- sum((all_draws - draw_struct_mean)^2, na.rm = TRUE)

      rhat_within <- NA_real_
      if (length(structs_with_multi_chains) > 0) {
        rw_vals <- vapply(structs_with_multi_chains, function(sk) {
          idx <- which(chain_primary_struct == sk & chain_stats$Unique_Structures == 1)
          calc_rhat(ch_vecs[idx])
        }, numeric(1))
        if (any(!is.na(rw_vals))) {
          rhat_within <- max(rw_vals, na.rm = TRUE)
        }
      }

      if (has_groups) {
        # 3-tier partitioning
        draw_group <- rep(chain_group_vec, each = n_iters)
        group_means <- tapply(all_draws, draw_group, mean, na.rm = TRUE)
        draw_group_mean <- group_means[as.character(draw_group)]

        ss_within_group  <- sum((draw_struct_mean - draw_group_mean)^2, na.rm = TRUE)
        ss_between_group <- sum((draw_group_mean - grand_mean)^2, na.rm = TRUE)

        if (ss_total > 1e-12) {
          pct_mcmc   <- max(0, min(100, round(ss_mcmc / ss_total * 100, 1)))
          pct_struct <- max(0, min(100, round(ss_within_group / ss_total * 100, 1)))
          pct_group  <- max(0, min(100, round(100 - pct_mcmc - pct_struct, 1)))
        } else {
          pct_mcmc <- 100.0; pct_struct <- 0.0; pct_group <- 0.0
        }

        rhat_within_grp <- list()
        for (g in unique_groups) {
          idx_g <- which(chain_group_vec == g)
          rhat_within_grp[[g]] <- if (length(idx_g) >= 2) calc_rhat(ch_vecs[idx_g]) else NA_real_
        }

        grp_draws <- lapply(unique_groups, function(g) unlist(ch_vecs[which(chain_group_vec == g)]))
        rhat_group <- calc_group_rhat(grp_draws)
        rhat_all   <- calc_rhat(ch_vecs)

        status_s_sens  <- paste0(s_term, "-Sensitive")
        status_g_sens  <- paste0(g_term, "-Sensitive")
        status_sg_sens <- paste0(s_term, "-Sensitive (within ", g_term, ")")

        status_str <- if (!is.na(rhat_within) && rhat_within > 1.05) {
          "Check MCMC"
        } else if (!is.na(rhat_group) && rhat_group > 1.05) {
          status_g_sens
        } else if (any(vapply(rhat_within_grp, function(v) !is.na(v) && v > 1.05, logical(1)))) {
          status_sg_sens
        } else if (!is.na(rhat_all) && rhat_all > 1.05) {
          status_s_sens
        } else {
          "Invariant"
        }

        row_data <- list(
          Parameter = p,
          Mean = round(mean(all_draws, na.rm = TRUE), 3),
          SD = round(stats::sd(all_draws, na.rm = TRUE), 3),
          CI_2.5 = round(unname(stats::quantile(all_draws, 0.025, na.rm = TRUE)), 3),
          CI_97.5 = round(unname(stats::quantile(all_draws, 0.975, na.rm = TRUE)), 3),
          Rhat_MCMC = round(rhat_within, 3)
        )
        for (g in unique_groups) {
          col_name <- paste0("Rhat_", gsub("[^A-Za-z0-9_]", "_", g))
          row_data[[col_name]] <- round(rhat_within_grp[[g]], 3)
        }
        row_data[[if (is_phylo) "Rhat_backbone" else "Rhat_group"]] <- round(rhat_group, 3)
        row_data[["Rhat_all"]] <- round(rhat_all, 3)
        row_data[["Pct_MCMC"]] <- pct_mcmc
        row_data[[paste0("Pct_", s_term)]] <- pct_struct
        row_data[[paste0("Pct_", g_term)]] <- pct_group
        row_data[["Status"]]   <- status_str
        as.data.frame(row_data, stringsAsFactors = FALSE)

      } else {
        # 2-tier partitioning
        ss_struct <- sum((draw_struct_mean - grand_mean)^2, na.rm = TRUE)
        if (ss_total > 1e-12) {
          pct_mcmc   <- max(0, min(100, round(ss_mcmc / ss_total * 100, 1)))
          pct_struct <- max(0, min(100, round(100 - pct_mcmc, 1)))
        } else {
          pct_mcmc <- 100.0; pct_struct <- 0.0
        }

        rhat_struct <- calc_rhat(ch_vecs)
        struct_means_val <- vapply(all_visited, function(sk) {
          idx <- which(chain_primary_struct == sk)
          if (length(idx) == 0) return(NA_real_)
          mean(unlist(ch_vecs[idx]), na.rm = TRUE)
        }, numeric(1))
        struct_sd <- if (sum(!is.na(struct_means_val)) > 1) stats::sd(struct_means_val, na.rm = TRUE) else 0.0

        status_sens <- paste0(s_term, "-Sensitive")
        status_str <- if (!is.na(rhat_within)) {
          if (rhat_within > 1.05) {
            "Check MCMC"
          } else if (!is.na(rhat_struct) && rhat_struct > 1.05) {
            status_sens
          } else {
            "Invariant"
          }
        } else {
          if (!is.na(rhat_struct) && rhat_struct > 1.05) status_sens else "Invariant"
        }

        rhat_col_name <- if (is_phylo) "Rhat_tree" else "Rhat_struct"
        pct_col_name  <- paste0("Pct_", s_term)
        sd_col_name   <- paste0(s_term, "_Mean_SD")

        row_df <- data.frame(
          Parameter = p,
          Mean = round(mean(all_draws, na.rm = TRUE), 3),
          SD = round(stats::sd(all_draws, na.rm = TRUE), 3),
          CI_2.5 = round(unname(stats::quantile(all_draws, 0.025, na.rm = TRUE)), 3),
          CI_97.5 = round(unname(stats::quantile(all_draws, 0.975, na.rm = TRUE)), 3),
          Rhat_within = round(rhat_within, 3),
          stringsAsFactors = FALSE
        )
        row_df[[rhat_col_name]] <- round(rhat_struct, 3)
        row_df[["Pct_MCMC"]]    <- pct_mcmc
        row_df[[pct_col_name]]  <- pct_struct
        row_df[[sd_col_name]]   <- round(struct_sd, 3)
        row_df[["Status"]]      <- status_str
        row_df
      }
    })
    param_diag_df <- do.call(rbind, diag_rows)
  }

  if (is_balanced_ensemble) {
    status <- "ENSEMBLE_BMA"
    summary_msg <- if (has_groups) {
      paste0(
        sprintf("Mode: Nested Multi-%s %s Ensemble (BMA across %d %ss: %s).\n",
                g_term, s_term, length(unique_groups), g_term_l, paste(unique_groups, collapse = ", ")),
        sprintf("%s Coverage: %d candidate %ss across %d %ss equally weighted (%.1f%% posterior draws per %s).\n",
                s_term, total_unique, s_term_l, length(unique_groups), g_term_l, 100 / total_unique, s_term_l),
        if (!is.null(param_diag_df) && any(!is.na(param_diag_df$Rhat_MCMC))) {
          sprintf("Max Within-%s Rhat (MCMC): %.3f | Max Between-%s Rhat: %.3f | Max Total Rhat: %.3f",
                  s_term, max(param_diag_df$Rhat_MCMC, na.rm = TRUE),
                  g_term, max(param_diag_df[[if (is_phylo) "Rhat_backbone" else "Rhat_group"]], na.rm = TRUE),
                  max(param_diag_df$Rhat_all, na.rm = TRUE))
        } else {
          sprintf("Tip: Use n.chains >= 2 to compute within-%s MCMC convergence alongside %s sensitivity.", s_term_l, g_term_l)
        }
      )
    } else {
      paste0(
        sprintf("Mode: Deterministic %s Ensemble (Bayesian Model Averaging).\n", s_term),
        sprintf("%s Coverage: %d %ss equally weighted (%.1f%% posterior draws per %s).\n",
                s_term, total_unique, s_term_l, 100 / total_unique, s_term_l),
        if (!is.null(param_diag_df) && any(!is.na(param_diag_df$Rhat_within))) {
          sprintf("Max Within-%s Rhat (MCMC convergence): %.3f | Max Between-%s Rhat (Sensitivity): %.3f",
                  s_term, max(param_diag_df$Rhat_within, na.rm = TRUE),
                  s_term, max(param_diag_df[[if (is_phylo) "Rhat_tree" else "Rhat_struct"]], na.rm = TRUE))
        } else {
          sprintf("Tip: Use n.chains >= 2 to compute within-%s Rhat alongside between-%s Rhat.", s_term_l, s_term_l)
        }
      )
    }
  } else if (total_jumps == 0) {
    status <- "LOCKED"
    n_tot_msg <- if (!is.null(model$n_trees_total)) sprintf(" (%d of %d candidate %ss visited)", total_unique, model$n_trees_total, s_term_l) else ""
    summary_msg <- paste0(
      sprintf("Mode: Discrete Categorical %s Sampling (method = 'dcat') — LOCKED.\n", s_term),
      sprintf("Within-chain %s jumping: 0 transitions observed (0.00%% jump rate)%s.\n", s_term_l, n_tot_msg),
      sprintf("Recommendation: Use structure_method = 'ensemble' (default in because >= 1.9.0)\n"),
      sprintf("to guarantee exact, balanced Bayesian Model Averaging across all candidate %ss.", s_term_l)
    )
  } else if (mean(chain_stats$Jump_Rate_Pct) >= 5) {
    status <- "ACTIVE_JUMPING"
    summary_msg <- paste0(
      sprintf("Active within-chain %s jumping observed (method = 'dcat')!\n", s_term_l),
      "Average jump rate: ", round(mean(chain_stats$Jump_Rate_Pct), 2), "%.\n",
      sprintf("Chains are actively exploring multiple %ss within single runs.", s_term_l)
    )
  } else {
    status <- "SLUGGISH_MIXING"
    summary_msg <- paste0(
      sprintf("Sluggish %s jumping observed (method = 'dcat').\n", s_term_l),
      "Average jump rate: ", round(mean(chain_stats$Jump_Rate_Pct), 2), "%.\n",
      sprintf("Consider using structure_method = 'ensemble' for balanced %s representation.", s_term_l)
    )
  }

  var_decomp_df <- NULL
  if (!is.null(param_diag_df) && "Pct_MCMC" %in% names(param_diag_df)) {
    var_cols <- if (has_groups) {
      c("Parameter", "SD", "Pct_MCMC", paste0("Pct_", s_term), paste0("Pct_", g_term))
    } else {
      c("Parameter", "SD", "Pct_MCMC", paste0("Pct_", s_term))
    }
    var_decomp_df <- param_diag_df[, var_cols, drop = FALSE]
    if (is_phylo) {
      var_decomp_df$Pct_Struct <- var_decomp_df[[paste0("Pct_", s_term)]]
      if (has_groups) var_decomp_df$Pct_Group <- var_decomp_df[[paste0("Pct_", g_term)]]
    }
  }

  if (!is.null(param_diag_df) && is_phylo) {
    param_diag_df$Pct_Struct <- param_diag_df$Pct_Tree
    param_diag_df$Rhat_struct <- param_diag_df$Rhat_tree
    if (has_groups) {
      param_diag_df$Pct_Group <- param_diag_df$Pct_Backbone
      param_diag_df$Rhat_group <- param_diag_df$Rhat_backbone
    }
  }

  out <- list(
    parameter = k_col,
    is_phylo = is_phylo,
    s_term = s_term,
    g_term = g_term,
    n_chains = n_chains,
    n_iters = n_iters,
    total_unique_structures = total_unique,
    visited_structures = all_visited,
    structure_weights_pct = struct_weights_pct,
    chain_stats = chain_stats,
    crosstab = crosstab_mat,
    param_diagnostics = param_diag_df,
    variance_decomposition = var_decomp_df,
    status = status,
    summary = summary_msg,
    chain_data = chain_data,
    group = group,
    unique_groups = unique_groups,
    has_groups = has_groups
  )
  class(out) <- "check_multicov"

  print(out)

  if (plot) {
    plot(out, ...)
  }

  invisible(out)
}

#' @export
print.check_multicov <- function(x, ...) {
  is_phylo <- isTRUE(x$is_phylo)
  s_term   <- x$s_term %||% (if (is_phylo) "Tree" else "Structure")
  g_term   <- x$g_term %||% (if (is_phylo) "Backbone" else "Group")

  banner_title <- if (is_phylo) "multiPhylo Convergence & Sensitivity Report" else "multiCov Convergence & Sensitivity Report"
  cat("\n====================================================================\n")
  cat(sprintf("%*s\n", floor((68 + nchar(banner_title))/2), banner_title))
  cat("====================================================================\n")
  cat(sprintf("%-30s: %s\n", paste(s_term, "Parameter Monitored"), x$parameter))
  cat(sprintf("%-30s: %d\n", "Number of Chains", x$n_chains))
  cat(sprintf("%-30s: %d\n", "Iterations per Chain", x$n_iters))
  cat(sprintf("%-30s: %d (IDs: %s)\n", paste("Total", paste0(s_term, "s"), "Represented"), x$total_unique_structures, paste(x$visited_structures, collapse = ", ")))
  if (isTRUE(x$has_groups)) {
    cat(sprintf("%-30s: %s\n", paste(paste0(g_term, "s"), "Represented"), paste(x$unique_groups, collapse = ", ")))
  }
  cat("\n")

  if (!identical(x$status, "ENSEMBLE_BMA") && !is.null(x$chain_stats)) {
    cat(sprintf("--- Per-Chain %s Exploration ---\n", s_term))
    print(x$chain_stats, row.names = FALSE)
    cat("\n")
  }

  cat(sprintf("--- Posterior Weight per %s (%%) ---\n", s_term))
  print(x$structure_weights_pct)

  if (!is.null(x$param_diagnostics) && nrow(x$param_diagnostics) > 0) {
    if (isTRUE(x$has_groups)) {
      cat("\n--- Parameter Diagnostics (Convergence & Sensitivity Rhat) ---\n")
    } else {
      cat(sprintf("\n--- Split Parameter Diagnostics (Within-%s MCMC vs. Between-%s Sensitivity) ---\n", s_term, s_term))
    }
    p_df <- x$param_diagnostics
    if (is_phylo) {
      p_df$Pct_Struct <- NULL
      p_df$Rhat_struct <- NULL
      p_df$Pct_Group <- NULL
      p_df$Rhat_group <- NULL
    }
    print(p_df, row.names = FALSE)
  }

  if (!is.null(x$variance_decomposition) && nrow(x$variance_decomposition) > 0) {
    cat("\n--- Uncertainty Variance Decomposition (% of Total Posterior Variance) ---\n")
    v_df <- x$variance_decomposition
    if (is_phylo) {
      v_df$Pct_Struct <- NULL
      v_df$Pct_Group <- NULL
    }
    print(v_df, row.names = FALSE)
    cat(sprintf("Note: Pct_MCMC is within-%s posterior variance (data-driven estimation uncertainty);\n", tolower(s_term)))
    cat(sprintf("      Pct_%s and Pct_%s reflect excess variance from phylogenetic uncertainty.\n", s_term, g_term))
  }

  cat("\n--- Diagnostic Assessment ---\n")
  cat(x$summary, "\n")
  cat("====================================================================\n\n")
  invisible(x)
}

#' @export
plot.check_multicov <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(1)
    suppressWarnings(par(old_par))
  }, add = TRUE)

  is_phylo  <- isTRUE(x$is_phylo)
  s_term    <- x$s_term %||% (if (is_phylo) "Tree" else "Structure")
  g_term    <- x$g_term %||% (if (is_phylo) "Backbone" else "Group")
  pct_s_col <- if ("Pct_Tree" %in% names(x$param_diagnostics)) "Pct_Tree" else "Pct_Struct"
  pct_g_col <- if ("Pct_Backbone" %in% names(x$param_diagnostics)) "Pct_Backbone" else "Pct_Group"
  rhat_s_col <- if ("Rhat_tree" %in% names(x$param_diagnostics)) "Rhat_tree" else "Rhat_struct"
  rhat_g_col <- if ("Rhat_backbone" %in% names(x$param_diagnostics)) "Rhat_backbone" else "Rhat_group"

  if (identical(x$status, "ENSEMBLE_BMA") && !is.null(x$param_diagnostics) && nrow(x$param_diagnostics) > 0) {
    param_diag_df <- x$param_diagnostics
    max_p_len <- max(nchar(param_diag_df$Parameter))
    left_mar <- max(5.5, min(9.0, max_p_len * 0.35 + 1.0))

    has_decomp <- all(c("Pct_MCMC", pct_s_col) %in% names(param_diag_df))
    if (has_decomp) {
      if (isTRUE(x$has_groups) && pct_g_col %in% names(param_diag_df)) {
        v_mat <- rbind(
          MCMC   = param_diag_df$Pct_MCMC,
          Struct = param_diag_df[[pct_s_col]],
          Group  = param_diag_df[[pct_g_col]]
        )
        v_labels <- c(paste0("Within-", s_term, " (Data)"), paste0("Between-", s_term, "s (within ", g_term, ")"), paste0("Between-", g_term, "s"))
        v_cols   <- c("#1B9E77", "#386CB0", "#D95F02")
      } else {
        v_mat <- rbind(
          MCMC   = param_diag_df$Pct_MCMC,
          Struct = param_diag_df[[pct_s_col]]
        )
        v_labels <- c(paste0("Within-", s_term, " (Data)"), paste0("Between-", s_term, "s"))
        v_cols   <- c("#1B9E77", "#D95F02")
      }
      colnames(v_mat) <- param_diag_df$Parameter
      bp <- as.numeric(barplot(v_mat, horiz = TRUE, plot = FALSE))
      y_idx <- bp
      y_lim <- c(0, max(bp) + 0.8)

      # Proportional 2-panel layout so left panel has room for parameter labels
      layout(matrix(1:2, nrow = 1), widths = c(1.2, 1.0))
      par(mar = c(4.2, left_mar, 4.8, 1.5), oma = c(0, 0, 2.0, 0))
    } else {
      v_mat <- NULL
      y_idx <- seq_len(nrow(param_diag_df))
      y_lim <- c(0.5, nrow(param_diag_df) + 0.5)
      par(mfrow = c(1, 1), mar = c(4.2, left_mar, 4.8, 2.5), oma = c(0, 0, 0, 0))
    }

    # Panel 1: Rhat Forest Plot
    if (isTRUE(x$has_groups)) {
      r_within <- param_diag_df$Rhat_MCMC
      r_grp    <- param_diag_df[[rhat_g_col]]
      grp_cols <- grep("^Rhat_", names(param_diag_df), value = TRUE)
      grp_cols <- setdiff(grp_cols, c("Rhat_MCMC", "Rhat_group", "Rhat_backbone", "Rhat_all"))

      all_r_vals <- unlist(param_diag_df[, c("Rhat_MCMC", rhat_g_col, "Rhat_all", grp_cols), drop = FALSE])
      x_max <- max(c(1.15, all_r_vals), na.rm = TRUE) * 1.05

      plot(r_grp, y_idx, xlim = c(0.98, x_max), ylim = y_lim,
           pch = 18, col = "#D95F02", cex = 1.4, yaxt = "n",
           xlab = expression(paste("Potential Scale Reduction Factor (", widehat(R), ")")),
           ylab = "", main = "", las = 1)
      mtext(if (has_decomp) "Convergence & Sensitivity (Rhat)" else paste0("Multi-", g_term, " Convergence vs. ", s_term, " Sensitivity"),
            side = 3, line = 3.2, font = 2, cex = 1.05)
      axis(2, at = y_idx, labels = param_diag_df$Parameter, las = 1, cex.axis = 0.75)
      grid()
      abline(v = 1.05, lty = 2, col = "firebrick", lwd = 1.5)

      if (any(!is.na(r_within))) {
        segments(r_within, y_idx, r_grp, y_idx, col = "gray75", lty = 1, lwd = 1.2)
        points(r_within, y_idx, pch = 19, col = "#1B9E77", cex = 1.2)
      }

      palette_cols <- c("#7570B3", "#386CB0", "#F0027F", "#BF5B17")
      for (g_i in seq_along(grp_cols)) {
        g_name <- grp_cols[g_i]
        col_val <- palette_cols[((g_i - 1) %% length(palette_cols)) + 1]
        points(param_diag_df[[g_name]], y_idx, pch = 17, col = col_val, cex = 1.2)
      }
      points(r_grp, y_idx, pch = 18, col = "#D95F02", cex = 1.4)

      leg_labels <- c(paste0("Within-", s_term, " (MCMC)"), sub("^Rhat_", "", grp_cols), paste0("Between-", g_term, "s"), "1.05 Cutoff")
      leg_cols   <- c("#1B9E77", palette_cols[seq_along(grp_cols)], "#D95F02", "firebrick")
      leg_pchs   <- c(19, rep(17, length(grp_cols)), 18, NA)
      leg_ltys   <- c(NA, rep(NA, length(grp_cols)), NA, 2)
      leg_lwds   <- c(NA, rep(NA, length(grp_cols)), NA, 1.5)

      legend("top", inset = c(0, -0.22), xpd = TRUE,
             ncol = if (length(grp_cols) > 1) 3 else length(leg_labels),
             legend = leg_labels, col = leg_cols, pch = leg_pchs,
             lty = leg_ltys, lwd = leg_lwds, pt.cex = 1.1, bty = "n", cex = 0.72)

    } else {
      r_within <- param_diag_df$Rhat_within
      r_struct <- param_diag_df[[rhat_s_col]]
      x_max <- max(c(1.15, r_within, r_struct), na.rm = TRUE) * 1.05

      plot(r_struct, y_idx, xlim = c(0.98, x_max), ylim = y_lim,
           pch = 17, col = "#D95F02", cex = 1.3, yaxt = "n",
           xlab = expression(paste("Potential Scale Reduction Factor (", widehat(R), ")")),
           ylab = "", main = "", las = 1)
      mtext(if (has_decomp) "Convergence & Sensitivity (Rhat)" else paste0("Ensemble Convergence vs. ", s_term, " Sensitivity"),
            side = 3, line = 3.2, font = 2, cex = 1.05)
      axis(2, at = y_idx, labels = param_diag_df$Parameter, las = 1, cex.axis = 0.75)
      grid()
      abline(v = 1.05, lty = 2, col = "firebrick", lwd = 1.5)
      if (any(!is.na(r_within))) {
        segments(r_within, y_idx, r_struct, y_idx, col = "gray70", lty = 1, lwd = 1.2)
        points(r_within, y_idx, pch = 19, col = "#1B9E77", cex = 1.2)
      }
      points(r_struct, y_idx, pch = 17, col = "#D95F02", cex = 1.3)

      legend("top", inset = c(0, -0.18), xpd = TRUE, horiz = TRUE,
             legend = c(paste0("Within-", s_term, " (MCMC)"), paste0("Between-", s_term, "s"), "1.05 Cutoff"),
             col = c("#1B9E77", "#D95F02", "firebrick"),
             pch = c(19, 17, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 1.5),
             pt.cex = 1.1, bty = "n", cex = 0.75)
    }

    # Panel 2: Stacked Variance Decomposition Bar Plot
    if (has_decomp && !is.null(v_mat)) {
      par(mar = c(4.2, 1.5, 4.8, 2.0))
      barplot(v_mat, horiz = TRUE, ylim = y_lim, col = v_cols,
              xlab = "Total Posterior Uncertainty (%)", xlim = c(0, 100), yaxt = "n",
              main = "", las = 1)
      mtext("Uncertainty Variance Decomposition (%)", side = 3, line = 3.2, font = 2, cex = 1.05)
      grid(nx = NULL, ny = NA)
      legend("top", inset = c(0, -0.18), xpd = TRUE, horiz = TRUE,
             legend = v_labels, fill = v_cols, bty = "n", cex = 0.78)

      banner_txt <- paste0(if (is_phylo) "multiPhylo" else "multiCov",
                           " Ensemble Diagnostics & Uncertainty Variance Decomposition")
      mtext(banner_txt, outer = TRUE, font = 2, cex = 1.1)
    }

  } else {
    # Categorical Jumping Diagnostic Panels
    layout(matrix(1:2, nrow = 1), widths = c(1.0, 1.0))
    par(mar = c(4.5, 4.5, 3.2, 1.2), oma = c(0, 0, 1.5, 0))
    base_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D")
    cols <- rep(base_palette, length.out = x$n_chains)

    # Panel 1: Trace plot of structure index over iterations
    y_range <- range(x$visited_structures)
    y_lim <- if (diff(y_range) == 0) {
      c(max(1, y_range[1] - 1), y_range[1] + 1)
    } else {
      c(max(1, y_range[1] - 0.5), y_range[2] + 0.5)
    }

    plot(1:x$n_iters, x$chain_data[[1]], type = "l", col = cols[1], ylim = y_lim,
         xlab = "Iteration", ylab = paste0(s_term, " Index (", x$parameter, ")"),
         main = paste("Trace of", s_term, "Selection"), lwd = 1.3, las = 1)
    grid()
    if (x$n_chains > 1) {
      for (ch in 2:x$n_chains) {
        lines(1:x$n_iters, x$chain_data[[ch]], col = cols[ch], lwd = 1.3)
      }
    }

    # Panel 2: Frequency barplot of visited structures by chain
    barplot(t(x$crosstab), beside = TRUE, col = cols,
            xlab = paste(s_term, "Index"), ylab = "Posterior Iterations",
            main = paste(s_term, "Frequencies per Chain"), las = 1)
    grid(nx = NA, ny = NULL)

    mtext(paste0(if (is_phylo) "multiPhylo" else "multiCov", " Diagnostic Report"),
          outer = TRUE, font = 2, cex = 1.1)
  }

  invisible(x)
}
