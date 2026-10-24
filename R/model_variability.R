#' Generate measurement error, variability, and predictor imputation code
#'
#' Adds observation models for variables with measurement error (se, reps),
#' and generates imputation priors for unmodelled predictors and latent variables.
#'
#' @param ctx Builder context environment
#' @keywords internal
generate_model_variability <- function(ctx) {
  # 1. Measurement error / Variability observation models
  if (length(ctx$variability_list) > 0) {
    ctx$model_lines <- safe_add_lines(ctx$model_lines, "  # Measurement error / Variability", ctx$declared_nodes)

    for (var in names(ctx$variability_list)) {
      type <- ctx$variability_list[[var]]

      var_dist <- ctx$dist_list[[var]] %||% "gaussian"
      if (!var_dist %in% c("gaussian", "normal")) {
        next
      }

      if (!var %in% ctx$all_vars) {
        warning(paste("Variable", var, "specified in 'variability' but not found in equations."))
        next
      }

      var_loop_N <- get_loop_bound(var, ctx$hierarchical_info, default_N = ctx$main_loop_N)

      if (type == "se") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("  for (i in 1:", var_loop_N, ") {"),
          paste0("    ", var, "_tau_obs[i] <- 1/(", var, "_se[i] * ", var, "_se[i])"),
          paste0("    ", var, "_mean[i] ~ dnorm(", var, "[i], ", var, "_tau_obs[i])")
        )
        if (ctx$engine == "jags") {
          ctx$model_lines <- c(
            ctx$model_lines,
            paste0("    log_lik_", var, "_mean[i] <- logdensity.norm(", var, "_mean[i], ", var, "[i], ", var, "_tau_obs[i])")
          )
        }
        ctx$model_lines <- c(ctx$model_lines, "  }")
      } else if (type == "reps") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("  for (i in 1:", var_loop_N, ") {"),
          paste0("    for (j in 1:N_reps_", var, "[i]) {"),
          paste0("      ", var, "_obs[i, j] ~ dnorm(", var, "[i], ", var, "_tau)")
        )
        if (ctx$engine == "jags") {
          ctx$model_lines <- c(
            ctx$model_lines,
            paste0("      lik_matrix_", var, "[i, j] <- logdensity.norm(", var, "_obs[i, j], ", var, "[i], ", var, "_tau)")
          )
        }
        ctx$model_lines <- c(ctx$model_lines, "    }")
        if (ctx$engine == "jags") {
          ctx$model_lines <- c(
            ctx$model_lines,
            paste0("    log_lik_", var, "_reps[i] <- sum(lik_matrix_", var, "[i, 1:N_reps_", var, "[i]])")
          )
        }
        ctx$model_lines <- c(ctx$model_lines, "  }")
        ctx$model_lines <- safe_add_lines(
          ctx$model_lines,
          c(
            get_precision_prior(paste0(var, "_tau"), var, priors = ctx$priors, family = ctx$family),
            paste0("  ", var, "_sigma <- 1/sqrt(", var, "_tau)")
          ),
          ctx$declared_nodes
        )
      } else {
        warning(paste("Unknown variability type:", type, "for variable", var))
      }
    }
  }

  # 2. Imputation priors for predictors (those not modeled as responses)
  non_response_vars <- setdiff(ctx$all_vars, names(ctx$response_counter))
  if (
    length(non_response_vars) > 0 &&
      (length(ctx$latent) > 0 ||
        (!is.null(ctx$vars_with_na) && length(ctx$vars_with_na) > 0))
  ) {
    ctx$model_lines <- safe_add_lines(ctx$model_lines, "  # Predictor priors for imputation", ctx$declared_nodes)
  }

  for (var in unique(non_response_vars)) {
    is_latent <- !is.null(ctx$latent) && var %in% ctx$latent

    if (!is_latent && (is.null(ctx$vars_with_na) || !var %in% ctx$vars_with_na)) {
      next
    }

    loop_bound_var <- get_loop_bound(var, ctx$hierarchical_info, default_N = ctx$main_loop_N)

    ctx$model_lines <- safe_add_lines(
      ctx$model_lines,
      c(
        paste0("  for (i in 1:", loop_bound_var, ") {"),
        paste0("    mu", var, "[i] <- 0"),
        paste0("  }")
      ),
      ctx$declared_nodes
    )

    if (ctx$independent) {
      if (is_latent && ctx$standardize_latent) {
        ctx$model_lines <- safe_add_lines(
          ctx$model_lines,
          c(
            paste0("  for (i in 1:", loop_bound_var, ") {"),
            paste0("    ", var, "[i] ~ dnorm(mu", var, "[i], 1)  # Standardized latent variable with parents"),
            paste0("  }")
          ),
          ctx$declared_nodes
        )
      } else {
        ctx$model_lines <- safe_add_lines(
          ctx$model_lines,
          c(
            paste0("  for (i in 1:", loop_bound_var, ") {"),
            paste0("    ", var, "[i] ~ dnorm(mu", var, "[i], tau_res_", var, ")"),
            paste0("  }")
          ),
          ctx$declared_nodes
        )
      }
    } else {
      if (is_latent && ctx$standardize_latent) {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("  for (i in 1:", loop_bound_var, ") {"),
          paste0("    ", var, "[i] ~ dnorm(mu", var, "[i], 1)  # Standardized latent variable with parents"),
          paste0("  }")
        )
      } else {
        additive_terms <- ""

        if (!is.null(ctx$structures)) {
          for (s_name in names(ctx$structures)) {
            if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(var, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) {
              next
            }
            
            s_suffix <- paste0("_", s_name)
            u_var <- paste0("u_", var, s_suffix)
            tau_u <- paste0("tau_u_", var, s_suffix)

            s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
            s_bound <- if (is.null(s_lvl)) loop_bound_var else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_idx_var <- get_struct_index(s_name, var, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

            r_lvl <- get_var_level(var, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
            local_struct_count <- 0
            for (sn in names(ctx$structures)) {
              if (is_valid_structure_mapping(get_struct_lvl(sn, ctx$hierarchical_info), r_lvl, ctx$hierarchical_info, allow_identity = TRUE)) {
                local_struct_count <- local_struct_count + 1
              }
            }
            use_partitioning <- (local_struct_count == 1) && 
                               (get_family_object(ctx$dist_list[[var]] %||% "gaussian")$family == "gaussian") &&
                               identical(s_lvl, r_lvl)

            s_def <- jags_structure_definition(
              ctx$structures[[s_name]],
              variable_name = var,
              s_name = s_name,
              loop_bound = s_bound,
              zeros_name = s_zeros,
              is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
              i_index = s_idx_var,
              use_partitioning = use_partitioning,
              engine = ctx$engine
            )

            if (!is.null(s_def$model_lines)) {
              ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
            }

            if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
               additive_terms <- paste0(additive_terms, " + ", s_def$term)
            }
          }
        }

        tau_e <- paste0("tau_res_", var)

        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("  for (i in 1:", loop_bound_var, ") {"),
          paste0(
            "    ",
            var,
            "[i] ~ dnorm(mu",
            var,
            "[i]",
            additive_terms,
            ", ",
            tau_e,
            ")"
          ),
          paste0("  }")
        )
      }
    }

    if (is_latent && !ctx$standardize_latent) {
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  # Latent variable: standardized (var = 1)"),
        paste0("  lambda", var, " ~ dunif(0, 1)"),
        paste0("  tau_u_", var, " <- 1/lambda", var),
        paste0("  tau_res_", var, " <- 1/(1-lambda", var, ")"),
        paste0("  sigma", var, " <- 1")
      )
    } else if (is_latent && ctx$standardize_latent) {
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  # Latent variable: fully standardized with N(0,1) prior (no variance parameters)")
      )
    } else {
      if (ctx$independent) {
        local_struct_count <- 0
        s_name_partition <- NULL
        if (!is.null(ctx$structures)) {
          for (sn in names(ctx$structures)) {
            if (is_valid_structure_mapping(get_struct_lvl(sn, ctx$hierarchical_info), get_var_level(var, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) {
              local_struct_count <- local_struct_count + 1
              s_name_partition <- sn
            }
          }
        }

        if (local_struct_count == 1 && !is.null(s_name_partition)) {
          partition_param <- paste0("lambda_", var)
          sigma_total_param <- paste0("sigma_total_", var)
          tau_total_param <- paste0("tau_total_", var)
          
          tau_struct_partition <- paste0("tau_u_", var, "_", s_name_partition)
          tau_res_partition <- paste0("tau_res_", var)

          ctx$model_lines <- safe_add_lines(
            ctx$model_lines,
            c(
              paste0("  # Variance Partitioning for independent var: ", var),
              paste0("  ", sigma_total_param, " ~ dunif(0, 10)"),
              paste0("  ", partition_param, " ~ dunif(0, 1)"),
              paste0("  ", tau_total_param, " <- 1/(", sigma_total_param, " * ", sigma_total_param, ")"),
              paste0("  ", tau_struct_partition, " <- ", tau_total_param, " / max(0.001, ", partition_param, ")"),
              paste0("  ", tau_res_partition, " <- ", tau_total_param, " / max(0.001, 1 - ", partition_param, ")")
            ),
            ctx$declared_nodes
          )
          ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = var, predictor = s_name_partition, parameter = partition_param, type = "structure")
          ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = var, predictor = s_name_partition, parameter = sigma_total_param, type = "structure")
        } else if (
          !is.null(ctx$fix_residual_variance) &&
            (var %in% names(ctx$fix_residual_variance) || length(ctx$fix_residual_variance) == 1)
        ) {
          val <- if (var %in% names(ctx$fix_residual_variance)) ctx$fix_residual_variance[[var]] else ctx$fix_residual_variance[[1]]
          prec <- 1 / val
          tau_line <- paste0("  tau_res_", var, " <- ", prec, " # Fixed residual variance")
          ctx$model_lines <- safe_add_lines(ctx$model_lines, tau_line, ctx$declared_nodes)
        } else {
          tau_line <- paste0("  ", get_precision_prior(paste0("tau_res_", var), var, priors = ctx$priors, family = ctx$family))
          ctx$model_lines <- safe_add_lines(ctx$model_lines, tau_line, ctx$declared_nodes)
        }
      } else {
        partition_already_handled <- FALSE
        s_name_handled <- NULL
        if (!is.null(ctx$structures)) {
          if (exists("s_def") && isTRUE(s_def$partition_handled) && s_def$variable_name == var) {
            partition_already_handled <- TRUE
            for (sn in names(ctx$structures)) {
              if (is_valid_structure_mapping(get_struct_lvl(sn, ctx$hierarchical_info), get_var_level(var, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) {
                s_name_handled <- sn
                break
              }
            }
          }
        }
        
        if (partition_already_handled) {
          ctx$param_map[[length(ctx$param_map) + 1]] <- list(
            response = var, 
            predictor = s_name_handled %||% "structure", 
            parameter = paste0("lambda_", var), 
            type = "structure"
          )
        } else {
          if (
            !is.null(ctx$fix_residual_variance) &&
              (var %in% names(ctx$fix_residual_variance) || length(ctx$fix_residual_variance) == 1)
          ) {
            val <- if (var %in% names(ctx$fix_residual_variance)) ctx$fix_residual_variance[[var]] else ctx$fix_residual_variance[[1]]
            prec <- 1 / val
            tau_line <- paste0("  tau_res_", var, " <- ", prec, " # Fixed residual variance")
          } else {
            tau_line <- paste0("  ", get_precision_prior(paste0("tau_res_", var), var, priors = ctx$priors, family = ctx$family))
          }
          ctx$model_lines <- safe_add_lines(ctx$model_lines, tau_line, ctx$declared_nodes)
        }
      }
    }
  }

  invisible(ctx)
}
