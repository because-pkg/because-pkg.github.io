#' Generate observation likelihoods for JAGS/BUGS model
#'
#' Constructs the observation likelihood statements for all response variables,
#' including structured errors (phylogenetic, spatial, multiPhylo),
#' unstructured random intercepts, and MAG correlated errors.
#'
#' Dispatches to per-family helpers defined in:
#'   - likelihood_continuous.R  (gaussian)
#'   - likelihood_categorical.R (multinomial, ordinal)
#'   - likelihood_discrete.R    (binomial, poisson, negbinomial/zinb, extensions)
#'
#' @param ctx Builder context environment
#' @keywords internal
generate_likelihoods <- function(ctx) {
  ctx$model_lines <- safe_add_lines(ctx$model_lines, "  # Multivariate normal likelihoods", ctx$declared_nodes)

  # Likelihoods for responses
  for (response in unique(names(ctx$response_counter))) {
    dist <- ctx$dist_list[[response]] %||% "gaussian"
    if (response %in% ctx$correlated_vars && dist == "gaussian") next

    for (k in 1:ctx$response_counter[[response]]) {
      suffix       <- if (k == 1) "" else as.character(k)
      dist         <- ctx$dist_list[[response]] %||% "gaussian"
      response_var <- paste0(response, suffix)
      loop_bound   <- get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N)
      eq_loop_N    <- loop_bound   # used by extension dispatch

      if (dist == "gaussian") {
        ctx <- jags_likelihood_gaussian(ctx, response, suffix, loop_bound)
      } else if (dist == "binomial") {
        ctx <- jags_likelihood_binomial(ctx, response, suffix, loop_bound)
      } else if (dist == "multinomial") {
        ctx <- jags_likelihood_multinomial(ctx, response, suffix, loop_bound)
      } else if (dist == "ordinal") {
        ctx <- jags_likelihood_ordinal(ctx, response, suffix, loop_bound)
      } else if (dist == "poisson") {
        ctx <- jags_likelihood_poisson(ctx, response, suffix, loop_bound)
      } else if (dist == "negbinomial" || dist == "zinb") {
        ctx <- jags_likelihood_negbinomial(ctx, response, suffix, loop_bound)
      } else {
        ctx <- jags_likelihood_extension(ctx, response, suffix, loop_bound, dist, eq_loop_N)
      }
    }
  }

  # Likelihoods for MAG correlated variables (gaussian only)
  for (var in names(ctx$vars_error_terms)) {
    dist <- if (!is.null(ctx$dist_list[[var]])) ctx$dist_list[[var]] else "gaussian"
    if (dist != "gaussian") next

    err_terms  <- ctx$vars_error_terms[[var]]
    suffix     <- if ((ctx$response_counter[[var]] %||% 0) > 1) "1" else ""
    loop_bound <- get_loop_bound(var, ctx$hierarchical_info, default_N = ctx$main_loop_N)

    structure_term_str <- ""
    structure_terms    <- character()
    if (length(ctx$structures) > 0) {
      for (s_idx in seq_along(ctx$structures)) {
        s_name <- names(ctx$structures)[s_idx]
        if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
        s_obj  <- ctx$structures[[s_idx]]

        local_struct_count <- 0
        for (sn in names(ctx$structures)) {
          if (is_valid_structure_mapping(
            get_struct_lvl(sn, ctx$hierarchical_info),
            get_var_level(var, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars),
            ctx$hierarchical_info, allow_identity = TRUE
          )) {
            local_struct_count <- local_struct_count + 1
          }
        }

        use_partitioning_flag <- (dist == "gaussian" && local_struct_count == 1 && length(err_terms) == 0)

        s_def <- jags_structure_definition(
          s_obj,
          variable_name    = var,
          s_name           = s_name,
          loop_bound       = loop_bound,
          is_multi         = is_struct_multi(s_name, ctx$hierarchical_info),
          use_partitioning = use_partitioning_flag,
          engine           = ctx$engine
        )

        if (!is.null(s_def$model_lines)) {
          ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
        }
        if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
          structure_terms <- c(structure_terms, s_def$term)
        }

        tau_u_var   <- paste0("tau_u_", var, "_", s_name)
        sigma_u_var <- paste0("sigma_", var, "_", s_name)
        ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = var, predictor = s_name, parameter = tau_u_var,   equation_index = NA, type = "structure")
        ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = var, predictor = s_name, parameter = sigma_u_var, equation_index = NA, type = "structure")
      }
      if (length(structure_terms) > 0) {
        structure_term_str <- paste0(" + ", paste(structure_terms, collapse = " + "))
      }
    }

    error_sum   <- paste(c(err_terms), collapse = " + ")
    mu_var      <- paste0("mu_", var, suffix)
    tau_res_var <- paste0("tau_res_", var, suffix)

    ctx$model_lines <- c(
      ctx$model_lines,
      paste0("  # Likelihood for ", var, " (with correlated residual errors)"),
      paste0("  for (i in 1:", loop_bound, ") {"),
      paste0("    ", var, "[i] ~ dnorm(", mu_var, "[i] + ", error_sum, structure_term_str, ", ", tau_res_var, ")")
    )
    if (ctx$engine == "jags") {
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    log_lik_", var, suffix, "[i] <- logdensity.norm(", var, "[i], ", mu_var, "[i] + ", error_sum, structure_term_str, ", ", tau_res_var, ")")
      )
    }
    ctx$model_lines <- safe_add_lines(ctx$model_lines, "  }", ctx$declared_nodes)
  }

  invisible(ctx)
}
