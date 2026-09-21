#' Generate preamble lines for JAGS/BUGS model
#'
#' Sets up zero vectors, Wishart priors for MAG induced correlations,
#' structural constants/priors, exogenous latent priors, and deterministic nodes.
#'
#' @param ctx Builder context environment
#' @keywords internal
generate_model_preamble <- function(ctx) {
  # Common structures and priors header
  ctx$model_lines <- c(
    "model {",
    "  # Common structures and priors"
  )

  # --- Handle Induced Correlations (MAG) Pair Processing ---
  if (!is.null(ctx$induced_correlations)) {
    ctx$model_lines <- c(
      ctx$model_lines,
      "  # Induced Correlations (Latent Variables) - Pair Processing"
    )

    for (pair in ctx$induced_correlations) {
      var1 <- pair[1]
      var2 <- pair[2]

      if (is.null(ctx$vars_error_terms[[var1]])) {
        ctx$vars_error_terms[[var1]] <- c()
      }
      if (is.null(ctx$vars_error_terms[[var2]])) {
        ctx$vars_error_terms[[var2]] <- c()
      }

      res_err <- paste0("err_res_", var1, "_", var2)
      tau_res_matrix <- paste0("TAU_res_", var1, "_", var2)
      cov_matrix <- paste0("cov_", var1, "_", var2)

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  # Correlated residuals between ", var1, " and ", var2, " (Wishart Prior)"),
        paste0("  ", tau_res_matrix, "[1:2, 1:2] ~ dwish(ID2[1:2, 1:2], 3)"),
        paste0("  ", cov_matrix, "[1:2, 1:2] <- inverse(", tau_res_matrix, "[1:2, 1:2])"),
        paste0("  sigma_res_", var1, "_", var2, " <- sqrt(", cov_matrix, "[1, 1])"),
        paste0("  sigma_res_", var2, "_", var1, " <- sqrt(", cov_matrix, "[2, 2])"),
        paste0("  rho_", var1, "_", var2, " <- ", cov_matrix, "[1, 2] / (sigma_res_", var1, "_", var2, " * sigma_res_", var2, "_", var1, ")")
      )

      loop_bound <- get_loop_bound(var1, ctx$hierarchical_info, default_N = ctx$main_loop_N)
      if (is.null(loop_bound) || nchar(trimws(loop_bound)) == 0) loop_bound <- "N"
      
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  for (i in 1:", loop_bound, ") {"),
        paste0("    ", res_err, "[i, 1:2] ~ dmnorm(zero_vec[1:2], ", tau_res_matrix, "[1:2, 1:2])"),
        paste0("  }")
      )

      ctx$vars_error_terms[[var1]] <- c(ctx$vars_error_terms[[var1]], paste0(res_err, "[i, 1]"))
      ctx$vars_error_terms[[var2]] <- c(ctx$vars_error_terms[[var2]], paste0(res_err, "[i, 2]"))
    }
  }

  ctx$model_lines <- safe_add_lines(ctx$model_lines, "  # Structural equations", ctx$declared_nodes)

  # --- Generic Structure Setup ---
  if (!is.null(ctx$structures)) {
    for (s_name in names(ctx$structures)) {
      def <- jags_structure_definition(
        ctx$structures[[s_name]],
        variable_name = "err",
        s_name = s_name,
        optimize = TRUE,
        engine = ctx$engine
      )
      if (!is.null(def$setup_code)) {
        ctx$model_lines <- safe_add_lines(ctx$model_lines, def$setup_code, ctx$declared_nodes)
      }
    }
  }

  # --- Handle Exogenous Latent Variables ---
  if (!is.null(ctx$variability)) {
    all_responses <- unique(sapply(ctx$equations, function(eq) {
      as.character(stats::formula(eq))[2]
    }))

    vars_with_variability <- names(ctx$variability_list)
    exogenous_vars <- setdiff(vars_with_variability, all_responses)

    if (ctx$is_multi_structure && is.null(ctx$structures)) {
      ctx$model_lines <- c(
        ctx$model_lines,
        "  # Multi-object sampling",
        "  K ~ dcat(p_obj[])",
        "  for (k in 1:Nobj) {",
        "    p_obj[k] <- 1/Nobj",
        "  }"
      )
    }

    if (length(exogenous_vars) > 0) {
      ctx$model_lines <- c(
        ctx$model_lines,
        "  # Priors for exogenous latent variables (variable with error but no parent)",
        "  for (i in 1:",
        if (!is.null(ctx$hierarchical_info)) ctx$hierarchical_info$counts[[1]] else "N",
        ") {"
      )
      for (ex_var in exogenous_vars) {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    ", ex_var, "[i] ~ dnorm(0, 1.0E-06)")
        )
      }
      ctx$model_lines <- c(ctx$model_lines, "  }")
    }
  }

  # --- Deterministic polynomial transformations ---
  if (!is.null(ctx$poly_terms) && length(ctx$poly_terms) > 0) {
    ctx$model_lines <- c(
      ctx$model_lines,
      "    # Deterministic polynomial transformations"
    )

    for (pt in ctx$poly_terms) {
      layout_N <- get_loop_bound(pt$base_var, ctx$hierarchical_info, default_N = ctx$main_loop_N)
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  for (i in 1:", layout_N, ") {"),
        sprintf("    %s[i] <- %s[i]^%d", pt$internal_name, pt$base_var, pt$power),
        "  }"
      )
    }
  }

  # --- Generic deterministic nodes (Interactions / Logic) ---
  if (!is.null(ctx$deterministic_terms) && length(ctx$deterministic_terms) > 0) {
    ctx$model_lines <- c(
      ctx$model_lines,
      "    # Deterministic nodes (Interactions / Logic)"
    )

    for (dt_name in names(ctx$deterministic_terms)) {
      dt <- ctx$deterministic_terms[[dt_name]]

      vars_in_expr_raw <- unique(gsub(
        "\\[i\\]",
        "",
        unlist(regmatches(
          dt$expression,
          gregexpr("\\b\\w+\\[i\\]", dt$expression)
        ))
      ))
      vars_in_expr <- unique(c(
        all.vars(parse(text = dt$original)),
        vars_in_expr_raw
      ))

      current_N <- ctx$main_loop_N
      finest_level_of_expr <- get_finest_level_of_vars(
        vars_in_expr,
        ctx$hierarchical_info,
        equations = ctx$equations,
        latent = ctx$latent,
        categorical_vars = ctx$categorical_vars
      )

      if (!is.null(finest_level_of_expr)) {
        current_N <- paste0("N_", finest_level_of_expr)
      }

      expr_hierarchical <- dt$expression
      if (!is.null(ctx$hierarchical_info) && !is.null(finest_level_of_expr)) {
        vars_sorted <- vars_in_expr[order(nchar(vars_in_expr), decreasing = TRUE)]
        for (v in vars_sorted) {
          correct_idx <- get_pred_index(
            v,
            finest_level_of_expr,
            ctx$hierarchical_info,
            equations = ctx$equations,
            latent = ctx$latent,
            categorical_vars = ctx$categorical_vars
          )
          expr_hierarchical <- gsub(
            paste0("\\b", v, "\\[i\\]"),
            correct_idx,
            expr_hierarchical
          )
        }
      }

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  for (i in 1:", current_N, ") {"),
        sprintf("    %s[i] <- %s", dt$internal_name, expr_hierarchical),
        "  }"
      )
    }
  }

  ctx$model_lines <- safe_add_lines(ctx$model_lines, "", ctx$declared_nodes)
  invisible(ctx)
}
