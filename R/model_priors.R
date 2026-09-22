#' Generate priors and variance partitioning for JAGS/BUGS model
#'
#' Constructs priors for intercepts (alpha), slopes (beta), residual precisions (tau),
#' structural components (tau_u, sigma), Pagel's lambda variance partitioning,
#' and family-specific auxiliary parameters.
#'
#' @param ctx Builder context environment
#' @keywords internal
generate_model_priors <- function(ctx) {
  ctx$model_lines <- safe_add_lines(ctx$model_lines, "  # Priors for structural parameters", ctx$declared_nodes)

  # Priors for alpha, lambda, tau, sigma
  for (response in unique(names(ctx$response_counter))) {
    dist <- ctx$dist_list[[response]] %||% "gaussian"
    if (response %in% ctx$correlated_vars && dist == "gaussian") {
      next
    }

    if (dist == "multinomial" || dist == "ordinal") {
      next
    }

    aux_flags <- list(skip_likelihood = FALSE, skip_variance = FALSE)
    if (grepl("^p_", response)) {
      target_dist <- ctx$dist_list[[sub("^p_", "", response)]] %||% "gaussian"
      fam_check_aux <- get_family_object(target_dist)
      aux_flags <- is_auxiliary_equation(fam_check_aux, response, names(ctx$dist_list))
    }

    for (k in 1:ctx$response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      alpha_name <- paste0("alpha_", response, suffix)
      non_identity_dists <- c(
        "binomial", "zip", "zinb", "bernoulli",
        "multinomial", "ordinal", "poisson", "negbinomial"
      )
      default_alpha <- "dnorm(0, 0.01)"
      is_non_identity <- dist %in% non_identity_dists || !(dist %in% c("gaussian", "normal"))
      if (is_non_identity || isTRUE(aux_flags$skip_likelihood)) {
        default_alpha <- "dnorm(0, 1)"
      }
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  ", get_prior(alpha_name, type = "alpha", default = default_alpha, priors = ctx$priors))
      )

      needs_residual_variance <- !is_non_identity && !isTRUE(aux_flags$skip_variance)
      needs_random_variance <- !ctx$independent &&
        (length(ctx$structure_names) > 0 || length(ctx$random_structure_names) > 0) &&
        !is_non_identity && !isTRUE(aux_flags$skip_variance)

      if (
        (is.null(ctx$vars_with_na) || !response %in% ctx$vars_with_na || TRUE) &&
          (needs_residual_variance || needs_random_variance) &&
          !isTRUE(aux_flags$skip_variance)
      ) {
        if (ctx$independent) {
          if (
            !is.null(ctx$fix_residual_variance) &&
              (response %in% names(ctx$fix_residual_variance) || length(ctx$fix_residual_variance) == 1)
          ) {
            val <- if (response %in% names(ctx$fix_residual_variance)) {
              ctx$fix_residual_variance[[response]]
            } else {
              ctx$fix_residual_variance[[1]]
            }
            prec <- 1 / val
            ctx$model_lines <- c(
              ctx$model_lines,
              paste0("  tau_res_", response, suffix, " <- ", prec, " # Fixed residual variance")
            )
          } else {
            ctx$model_lines <- safe_add_lines(
              ctx$model_lines,
              paste0("  ", get_precision_prior(paste0("tau_res_", response, suffix), response, priors = ctx$priors, family = ctx$family)),
              ctx$declared_nodes
            )
          }
        } else {
          if (
            !is.null(ctx$fix_residual_variance) &&
              (response %in% names(ctx$fix_residual_variance) || length(ctx$fix_residual_variance) == 1)
          ) {
            val <- if (response %in% names(ctx$fix_residual_variance)) {
              ctx$fix_residual_variance[[response]]
            } else {
              ctx$fix_residual_variance[[1]]
            }
            prec <- 1 / val
            ctx$model_lines <- c(
              ctx$model_lines,
              paste0("  tau_res_", response, suffix, " <- ", prec, " # Fixed residual variance")
            )
          } else {
            if (!dist %in% c("negbinomial", "zinb")) {
              ctx$model_lines <- safe_add_lines(
                ctx$model_lines,
                paste0("  ", get_precision_prior(paste0("tau_res_", response, suffix), response, priors = ctx$priors, family = ctx$family)),
                ctx$declared_nodes
              )
            }
          }

          processed_struct_signals <- character(0)
          for (s_name in ctx$structure_names) {
            if (s_name %in% processed_struct_signals) next
            processed_struct_signals <- c(processed_struct_signals, s_name)

            if (
              !is_valid_structure_mapping(
                get_struct_lvl(s_name, ctx$hierarchical_info),
                get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars),
                ctx$hierarchical_info,
                allow_identity = TRUE
              )
            ) {
              next
            }
            s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
            s_suffix <- paste0("_", s_name)

            is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
            if (is_unified) {
              tau_u <- paste0("tau_u_", s_name, "_", response)
              sigma_name <- paste0("sigma_", s_name, "_", response)
            } else {
              tau_u <- paste0("tau_u_", response, suffix, s_suffix)
              sigma_name <- paste0("sigma_", response, suffix, s_suffix)
            }

            ctx$model_lines <- safe_add_lines(
              ctx$model_lines,
              c(
                paste0("  ", get_precision_prior(tau_u, response, priors = ctx$priors, family = ctx$family)),
                paste0("  ", sigma_name, " <- 1/sqrt(", tau_u, ")")
              ),
              ctx$declared_nodes
            )
          }

          if (
            length(ctx$structure_names) == 1 &&
              length(ctx$random_structure_names) == 0 &&
              !dist %in% c("negbinomial", "zinb") &&
              !response %in% ctx$partitioned_responses
          ) {
            s_name <- ctx$structure_names[1]
            s_suffix <- paste0("_", s_name)
            is_unified_l <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
            
            if (is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), response, ctx$hierarchical_info, allow_identity = TRUE)) {
              if (is_unified_l) {
                tau_u_name <- paste0("tau_u_", s_name, "_", response)
              } else {
                tau_u_name <- paste0("tau_u_", response, suffix, s_suffix)
              }

              ctx$model_lines <- safe_add_lines(
                ctx$model_lines,
                paste0(
                  "  lambda", response, suffix,
                  " <- (1/", tau_u_name, ") / ((1/", tau_u_name, ") + (1/tau_res_", response, suffix, "))"
                ),
                ctx$declared_nodes
              )
            }
          }

          processed_rand_signals <- character(0)
          for (r_name in ctx$random_structure_names) {
            if (r_name %in% processed_rand_signals) next
            processed_rand_signals <- c(processed_rand_signals, r_name)

            is_requested <- FALSE
            for (rt in ctx$random_terms) {
              if (identical(rt$response, response) && identical(rt$group, r_name)) {
                is_requested <- TRUE
                break
              }
            }
            if (!is_requested) next

            if (!is_valid_random_level(response, r_name, ctx$hierarchical_info, family = ctx$family, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars, quiet = ctx$quiet)) next
            
            s_suffix <- paste0("_", r_name)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            ctx$model_lines <- safe_add_lines(
              ctx$model_lines,
              c(
                paste0("  ", get_precision_prior(tau_u, response, priors = ctx$priors, family = ctx$family)),
                paste0("  sigma_", response, suffix, "_", r_name, " <- 1/sqrt(", tau_u, ")")
              ),
              ctx$declared_nodes
            )
          }
        }
      }
    }
  }

  # Priors for multinomial parameters (arrays)
  for (response in unique(names(ctx$response_counter))) {
    dist <- ctx$dist_list[[response]] %||% "gaussian"
    if (dist == "multinomial") {
      K_var <- paste0("K_", response)
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  # Priors for ", response, " (Multinomial)"),
        paste0("  alpha_", response, "[1] <- 0"),
        paste0("  for (k in 2:", K_var, ") {"),
        paste0("    alpha_", response, "[k] ~ dnorm(0, 0.01)"),
        "  }"
      )

      if (!ctx$independent) {
        if (length(ctx$structure_names) > 0) {
          for (s_name in ctx$structure_names) {
            if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next
            
            is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
            if (is_unified) {
              tau_u <- paste0("tau_u_", s_name, "_", response)
              sigma_name <- paste0("sigma_", s_name, "_", response)
            } else {
              tau_u <- paste0("tau_u_", response, "_", s_name)
              sigma_name <- paste0("sigma_", response, "_", s_name)
            }

            ctx$model_lines <- c(
              ctx$model_lines,
              paste0("  for (k in 2:", K_var, ") {"),
              paste0("    ", tau_u, "[k] ~ dgamma(0.001, 0.001)"),
              paste0("    ", sigma_name, "[k] <- 1/sqrt(", tau_u, "[k])"),
              "  }"
            )
          }
        }

        for (r_name in ctx$random_structure_names) {
          is_requested <- FALSE
          for (rt in ctx$random_terms) {
            if (identical(rt$response, response) && identical(rt$group, r_name)) {
              is_requested <- TRUE
              break
            }
          }
          if (!is_requested) next
          if (!is_valid_random_level(response, r_name, ctx$hierarchical_info, family = ctx$family, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars, quiet = ctx$quiet)) next

          s_suffix <- paste0("_", r_name)
          tau_u <- paste0("tau_u_", response, s_suffix)
          sigma_name <- paste0("sigma_", response, "_", r_name)

          ctx$model_lines <- c(
            ctx$model_lines,
            paste0("  for (k in 2:", K_var, ") {"),
            paste0("    ", tau_u, "[k] ~ dgamma(0.001, 0.001)"),
            paste0("    ", sigma_name, "[k] <- 1/sqrt(", tau_u, "[k])"),
            "  }"
          )
        }
      }

      for (eq in ctx$eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            beta_name <- paste0("beta_", response, "_", pred)
            ctx$model_lines <- c(
              ctx$model_lines,
              paste0("  ", beta_name, "[1] <- 0"),
              paste0("  for (k in 2:", K_var, ") {"),
              paste0("    ", get_prior(paste0(beta_name, "[k]"), type = "beta", priors = ctx$priors)),
              "  }"
            )
          }
        }
      }
    }
  }

  # Priors for ordinal parameters (cutpoints + variance components)
  for (response in unique(names(ctx$response_counter))) {
    dist <- ctx$dist_list[[response]] %||% "gaussian"
    if (dist == "ordinal") {
      K_var <- paste0("K_", response)

      for (k in 1:ctx$response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)

        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("  # Priors for ", response, suffix, " (Ordinal)"),
          paste0("  cutpoint_raw_", response, suffix, "[1] ~ dnorm(0, 0.1)"),
          paste0("  cutpoint_", response, suffix, "[1] <- cutpoint_raw_", response, suffix, "[1]"),
          paste0("  for (k in 2:(", K_var, "-1)) {"),
          paste0("    cutpoint_raw_", response, suffix, "[k] ~ dnorm(0, 0.1)"),
          paste0("    cutpoint_", response, suffix, "[k] <- cutpoint_", response, suffix, "[k-1] + exp(cutpoint_raw_", response, suffix, "[k])"),
          "  }"
        )

        if (!ctx$independent) {
          if (length(ctx$structure_names) > 0) {
            for (s_name in ctx$structure_names) {
              if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next

              is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
              if (is_unified) {
                tau_u <- paste0("tau_u_", s_name, "_", response)
                sigma_name <- paste0("sigma_", s_name, "_", response)
              } else {
                tau_u <- paste0("tau_u_", response, "_", s_name)
                sigma_name <- paste0("sigma_", response, "_", s_name)
              }

              ctx$model_lines <- c(
                ctx$model_lines,
                paste0("  ", tau_u, " ~ dgamma(0.001, 0.001)"),
                paste0("  ", sigma_name, " <- 1/sqrt(", tau_u, ")")
              )
            }
          }

          for (r_name in ctx$random_structure_names) {
            is_requested <- FALSE
            for (rt in ctx$random_terms) {
              if (identical(rt$response, response) && identical(rt$group, r_name)) {
                is_requested <- TRUE
                break
              }
            }
            if (!is_requested) next
            if (!is_valid_random_level(response, r_name, ctx$hierarchical_info, family = ctx$family, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars, quiet = ctx$quiet)) next

            s_suffix <- paste0("_", r_name)
            tau_u <- paste0("tau_u_", response, s_suffix)
            sigma_name <- paste0("sigma_", response, "_", r_name)

            ctx$model_lines <- c(
              ctx$model_lines,
              paste0("  ", tau_u, " ~ dgamma(0.001, 0.001)"),
              paste0("  ", sigma_name, " <- 1/sqrt(", tau_u, ")")
            )
          }
        }
      }

      for (eq in ctx$eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            beta_name <- paste0("beta_", response, "_", pred)
            ctx$model_lines <- c(
              ctx$model_lines,
              paste0("  ", get_prior(beta_name, type = "beta", priors = ctx$priors))
            )
          }
        }
      }
    }
  }

  # Priors for Negative Binomial parameters (size parameter r)
  for (response in unique(names(ctx$response_counter))) {
    dist <- ctx$dist_list[[response]] %||% "gaussian"
    if (dist %in% c("negbinomial", "zinb")) {
      for (k in 1:ctx$response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)
        r_name <- paste0("r_", response, suffix)
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("  # Priors for ", response, suffix, " (Negative Binomial size parameter)"),
          paste0("  ", r_name, " ~ dgamma(0.01, 0.01)")
        )
        ctx$param_map[[length(ctx$param_map) + 1]] <- list(
          response = response,
          predictor = "(Dispersion)",
          parameter = r_name,
          equation_index = NA,
          type = "structure"
        )
      }
    }
  }

  # Priors for correlated vars alphas (intercepts - handled here for Gaussian only)
  for (var in ctx$correlated_vars) {
    dist <- if (!is.null(ctx$dist_list[[var]])) ctx$dist_list[[var]] else "gaussian"
    if (dist == "gaussian") {
      alpha_var <- paste0("alpha_", var)
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("  ", get_prior(alpha_var, type = "alpha", priors = ctx$priors))
      )

      if (ctx$independent) {
        ctx$model_lines <- safe_add_lines(
          ctx$model_lines,
          paste0("  ", get_precision_prior(paste0("tau_res_", var), var, priors = ctx$priors, family = ctx$family)),
          ctx$declared_nodes
        )
      } else {
        if (length(ctx$structures) > 0) {
          for (s_idx in seq_along(ctx$structures)) {
            s_name <- names(ctx$structures)[s_idx]
            if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
            
            if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(var, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next
            
            tau_u_var <- paste0("tau_u_", var, "_", s_name)
            ctx$model_lines <- safe_add_lines(
              ctx$model_lines,
              paste0("  ", get_precision_prior(tau_u_var, var, priors = ctx$priors, family = ctx$family)),
              ctx$declared_nodes
            )
          }
        }
        ctx$model_lines <- safe_add_lines(
          ctx$model_lines,
          paste0("  ", get_precision_prior(paste0("tau_res_", var), var, priors = ctx$priors, family = ctx$family)),
          ctx$declared_nodes
        )
      }
    }
  }

  # Priors for Zero-Inflation parameters (psi)
  for (response in unique(names(ctx$response_counter))) {
    dist <- ctx$dist_list[[response]] %||% "gaussian"
    if (dist %in% c("zip", "zinb")) {
      for (k in 1:ctx$response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)
        psi_name <- paste0("psi_", response, suffix)
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("  ", psi_name, " ~ dunif(0, 1) # Zero-inflation probability")
        )
      }
    }
  }

  # Priors for regression coefficients (beta)
  unique_betas <- unique(unlist(ctx$beta_counter))
  excluded_betas <- c()

  for (response in unique(names(ctx$response_counter))) {
    dist <- ctx$dist_list[[response]] %||% "gaussian"
    if (grepl("^p_", response)) {
      target_resp <- sub("^p_", "", response)
      target_dist <- ctx$dist_list[[target_resp]] %||% "gaussian"
      fam_check_aux <- get_family_object(target_dist)
      aux_flags <- is_auxiliary_equation(fam_check_aux, response, names(ctx$dist_list))
      if (isTRUE(aux_flags$skip_likelihood)) next
    }

    if (dist %in% c("multinomial", "ordinal")) {
      for (eq in ctx$eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            excluded_betas <- c(excluded_betas, paste0("beta_", response, "_", pred))
          }
        }
      }
    }
  }

  unique_betas <- setdiff(unique_betas, excluded_betas)
  pinned_latents <- character()

  for (beta in unique_betas) {
    default_beta <- "dnorm(0, 0.01)"

    found_resp <- NULL
    for (nm in names(ctx$dist_list)) {
      if (startsWith(beta, paste0("beta_", nm, "_"))) {
        found_resp <- nm
        break
      }
      if (startsWith(beta, paste0("beta_p_", nm, "_"))) {
        default_beta <- "dnorm(0, 1)"
        break
      }
    }

    if (!is.null(found_resp)) {
      d <- ctx$dist_list[[found_resp]] %||% "gaussian"
      if (d %in% non_identity_dists || !(d %in% c("gaussian", "normal"))) {
        default_beta <- "dnorm(0, 1)"
      }
    }

    if (!is.null(ctx$latent) && length(ctx$latent) > 0) {
      for (lat in ctx$latent) {
        if (!lat %in% pinned_latents && grepl(paste0("_", lat, "$"), beta)) {
          if (ctx$fix_latent == "loading") {
            default_beta <- "1.0"
          } else if (ctx$fix_latent == "sign") {
            default_beta <- paste0(default_beta, " T(0,)")
          }
          pinned_latents <- c(pinned_latents, lat)
          break
        }
      }
    }

    ctx$model_lines <- c(
      ctx$model_lines,
      paste0("  ", get_prior(beta, type = "beta", default = default_beta, priors = ctx$priors))
    )
  }

  if (!is.null(ctx$priors)) {
    community_extras <- ctx$priors[startsWith(names(ctx$priors), "__community__")]
    if (length(community_extras) > 0) {
      ctx$model_lines <- c(
        ctx$model_lines,
        "  # Community hyperparameter declarations (from extension package)"
      )
      for (line in unlist(community_extras)) {
        ctx$model_lines <- safe_add_lines(ctx$model_lines, paste0("  ", line), ctx$declared_nodes)
      }
    }
  }

  if (ctx$is_multi_structure) {
    ctx$model_lines <- c(
      ctx$model_lines,
      "  # Phylogenetic uncertainty weighting",
      "  for (k in 1:Ntree) {",
      "    p_tree[k] <- 1/Ntree",
      "  }",
      "  K ~ dcat(p_tree[1:Ntree])"
    )
  }

  invisible(ctx)
}
