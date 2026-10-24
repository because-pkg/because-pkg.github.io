# Categorical Likelihood Blocks for JAGS/NIMBLE Model Compiler
#
# Generates observation likelihood blocks for unordered (multinomial) and
# ordered (ordinal) categorical response variables, including random intercepts
# and phylogenetic/spatial structure terms.
#
# Primary functions:
#   - jags_likelihood_multinomial()
#   - jags_likelihood_ordinal()
#
# @keywords internal

jags_likelihood_multinomial <- function(ctx, response, suffix, loop_bound) {
  K_var <- paste0("K_", response)
  err   <- paste0("err_", response, suffix)

  if (ctx$independent) {
    ctx$model_lines <- c(
      ctx$model_lines,
      paste0("  # Independent Multinomial (no correlated error terms)"),
      paste0("  for (k in 2:", K_var, ") {"),
      paste0("    for (i in 1:", loop_bound, ") {"),
      paste0("      ", err, "[i, k] <- 0"),
      "    }",
      "  }"
    )
  } else {
    ctx$model_lines <- c(
      ctx$model_lines,
      paste0("  # Random effects for Multinomial: ", response),
      paste0("  for (k in 2:", K_var, ") {")
    )

    total_u <- ""
    if (length(ctx$structures) > 0) {
      for (s_idx in seq_along(ctx$structures)) {
        s_name <- names(ctx$structures)[s_idx]
        if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
        s_obj  <- ctx$structures[[s_idx]]

        if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next

        s_lvl     <- get_struct_lvl(s_name, ctx$hierarchical_info)
        s_bound   <- if (is.null(s_lvl)) get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N) else paste0("N_", s_lvl)
        s_zeros   <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
        s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

        s_def <- jags_structure_definition(
          s_obj,
          variable_name  = response,
          s_name         = s_name,
          loop_bound     = s_bound,
          zeros_name     = s_zeros,
          category_index = "k",
          is_multi       = is_struct_multi(s_name, ctx$hierarchical_info),
          i_index        = s_idx_var,
          engine         = ctx$engine
        )

        if (!is.null(s_def)) {
          ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
          total_u         <- paste0(total_u, " + ", s_def$term)

          is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
          if (is_unified) {
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", s_name, "_", response, "[k]"),   equation_index = NA, type = "structure")
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", s_name, "_", response, "[k]"), equation_index = NA, type = "structure")
          } else {
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name, "[k]"),   equation_index = NA, type = "structure")
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name, "[k]"), equation_index = NA, type = "structure")
          }
        }
      }
    }

    for (r_name in ctx$random_structure_names) {
      is_requested <- FALSE
      for (rt in ctx$random_terms) {
        if (identical(rt$response, response) && identical(rt$group, r_name)) { is_requested <- TRUE; break }
      }
      if (!is_requested) next
      if (!is_valid_random_level(response, r_name, ctx$hierarchical_info, family = ctx$family, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars, quiet = ctx$quiet)) next
      if (r_name %in% names(ctx$structures) || any(vapply(names(ctx$structures), function(s) grepl(s, r_name), logical(1)))) {
        is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
        if (is_unified_r) next
      }

      s_suffix_r <- paste0("_", r_name)
      u_std      <- paste0("u_std_", response, suffix, s_suffix_r)
      u          <- paste0("u_", response, suffix, s_suffix_r)
      tau_u      <- paste0("tau_u_", response, suffix, s_suffix_r)
      n_groups   <- paste0("N_", r_name)
      zeros_name <- paste0("zeros_", r_name)
      prec_name  <- paste0("Prec_", r_name)

      ctx$model_lines <- safe_add_lines(ctx$model_lines, c(
        paste0("  ", u_std, "[1:", n_groups, ", k] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
        paste0("  for (g in 1:", n_groups, ") { ", u, "[g, k] <- (", u_std, "[g, k] - mean(", u_std, "[1:", n_groups, ", k])) / sqrt(", tau_u, "[k]) }")
      ), ctx$declared_nodes)

      group_idx <- get_group_idx_string(response, r_name, ctx$hierarchical_info, default_N = ctx$main_loop_N)
      total_u   <- paste0(total_u, " + ", u, "[", group_idx, ", k]")

      ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0(tau_u, "[k]"),                             equation_index = NA, type = "structure")
      ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name, "[k]"), equation_index = NA, type = "structure")
    }

    ctx$model_lines <- c(
      ctx$model_lines,
      paste0("    for (i in 1:", loop_bound, ") {"),
      paste0("      ", err, "[i, k] <- ", if (nchar(total_u) > 0) sub("^ \\+ ", "", total_u) else "0"),
      "    }",
      "  }"
    )
  }

  invisible(ctx)
}


jags_likelihood_ordinal <- function(ctx, response, suffix, loop_bound) {
  err <- paste0("err_", response, suffix)

  if (ctx$independent) {
    ctx$model_lines <- c(
      ctx$model_lines,
      paste0("  # Independent Ordinal (no correlated error terms)"),
      paste0("  for (i in 1:", loop_bound, ") {"),
      paste0("    ", err, "[i] <- 0"),
      "  }"
    )
  } else {
    total_u <- ""

    if (length(ctx$structures) > 0) {
      for (s_idx in seq_along(ctx$structures)) {
        s_name <- names(ctx$structures)[s_idx]
        if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
        s_obj  <- ctx$structures[[s_idx]]

        if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next

        s_lvl     <- get_struct_lvl(s_name, ctx$hierarchical_info)
        s_bound   <- if (is.null(s_lvl)) get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N) else paste0("N_", s_lvl)
        s_zeros   <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
        s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

        s_def <- jags_structure_definition(
          s_obj,
          variable_name = response,
          s_name        = s_name,
          loop_bound    = s_bound,
          zeros_name    = s_zeros,
          is_multi      = is_struct_multi(s_name, ctx$hierarchical_info),
          i_index       = s_idx_var,
          engine        = ctx$engine
        )

        if (!is.null(s_def)) {
          ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
          total_u         <- paste0(total_u, " + ", s_def$term)

          is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
          if (is_unified) {
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", s_name, "_", response),   equation_index = NA, type = "structure")
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", s_name, "_", response), equation_index = NA, type = "structure")
          } else {
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name),   equation_index = NA, type = "structure")
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
          }
        }
      }
    }

    for (r_name in ctx$random_structure_names) {
      is_requested <- FALSE
      for (rt in ctx$random_terms) {
        if (identical(rt$response, response) && identical(rt$group, r_name)) { is_requested <- TRUE; break }
      }
      if (!is_requested) next
      if (!is_valid_random_level(response, r_name, ctx$hierarchical_info, family = ctx$family, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars, quiet = ctx$quiet)) next
      if (r_name %in% names(ctx$structures) || any(vapply(names(ctx$structures), function(s) grepl(s, r_name), logical(1)))) {
        is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
        if (is_unified_r) next
      }

      s_suffix_r <- paste0("_", r_name)
      u_std      <- paste0("u_std_", response, suffix, s_suffix_r)
      u          <- paste0("u_", response, suffix, s_suffix_r)
      tau_u      <- paste0("tau_u_", response, suffix, s_suffix_r)
      n_groups   <- paste0("N_", r_name)
      zeros_name <- paste0("zeros_", r_name)
      prec_name  <- paste0("Prec_", r_name)

      ctx$model_lines <- safe_add_lines(ctx$model_lines, c(
        paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
        paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- (", u_std, "[g] - mean(", u_std, "[1:", n_groups, "])) / sqrt(", tau_u, ") }")
      ), ctx$declared_nodes)

      group_idx <- get_group_idx_string(response, r_name, ctx$hierarchical_info, default_N = ctx$main_loop_N)
      total_u   <- paste0(total_u, " + ", u, "[", group_idx, "]")

      ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = r_name, parameter = tau_u, equation_index = NA, type = "structure")
      ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name), equation_index = NA, type = "structure")
    }

    mag_terms <- ctx$vars_error_terms[[response]]
    if (!is.null(mag_terms) && length(mag_terms) > 0) {
      total_u <- paste0(total_u, " + ", paste(mag_terms, collapse = " + "))
    }

    ctx$model_lines <- c(
      ctx$model_lines,
      paste0("  for (i in 1:", loop_bound, ") {"),
      paste0("    ", err, "[i] <- ", if (nchar(total_u) > 0) sub("^ \\+ ", "", total_u) else "0"),
      "  }"
    )
  }

  invisible(ctx)
}
