#' Generate observation likelihoods for JAGS/BUGS model
#'
#' Constructs the observation likelihood statements for all response variables,
#' including structured errors (phylogenetic, spatial, multiPhylo),
#' unstructured random intercepts, and MAG correlated errors.
#'
#' @param ctx Builder context environment
#' @keywords internal
generate_likelihoods <- function(ctx) {
  ctx$model_lines <- safe_add_lines(ctx$model_lines, "  # Multivariate normal likelihoods", ctx$declared_nodes)

  # Likelihoods for responses
  for (response in unique(names(ctx$response_counter))) {
    dist <- ctx$dist_list[[response]] %||% "gaussian"
    if (response %in% ctx$correlated_vars && dist == "gaussian") {
      next
    }

    dist <- ctx$dist_list[[response]] %||% "gaussian"

    for (k in 1:ctx$response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      dist <- ctx$dist_list[[response]] %||% "gaussian"
      response_var <- paste0(response, suffix)

      loop_bound <- get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N)
      resp_data_idx <- paste0(response, "[i]")
      tau <- paste0("TAU", tolower(response), suffix)

      if (dist == "gaussian") {
        if (grepl("^p_", response)) {
          all_resp_names <- names(ctx$dist_list)
          fam_check <- get_family_object(ctx$dist_list[[sub("^p_", "", response)]] %||% "gaussian")
          aux_flags <- is_auxiliary_equation(fam_check, response, all_resp_names)
          if (isTRUE(aux_flags$skip_likelihood)) next
        }

        mu <- paste0("mu_", response, suffix)
        tau_scalar <- paste0("tau", response, suffix)

        if (ctx$independent) {
          tau_res <- paste0("tau_res_", response, suffix)

          ctx$model_lines <- c(
            ctx$model_lines,
            paste0("  for (i in 1:", loop_bound, ") {"),
            paste0("    ", response_var, "[i] ~ dnorm(", mu, "[i], ", tau_res, ")")
          )
          if (ctx$engine == "jags") {
            ctx$model_lines <- c(
              ctx$model_lines,
              paste0("    log_lik_", response, suffix, "[i] <- logdensity.norm(", response_var, "[i], ", mu, "[i], ", tau_res, ")")
            )
          }
          ctx$model_lines <- safe_add_lines(ctx$model_lines, "  }", ctx$declared_nodes)
        } else {
          additive_terms <- ""
          processed_signals <- character(0)
          for (s_name in ctx$structure_names) {
            if (s_name %in% processed_signals) next
            processed_signals <- c(processed_signals, s_name)
            
            if (!is_valid_structure_mapping(
              get_struct_lvl(s_name, ctx$hierarchical_info),
              get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars),
              ctx$hierarchical_info,
              allow_identity = TRUE
            )) {
              next
            }

            s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
            s_bound <- if (is.null(s_lvl)) loop_bound else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

            r_lvl <- get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
            
            local_struct_count <- 0
            for (sn in ctx$structure_names) {
              if (is_valid_structure_mapping(get_struct_lvl(sn, ctx$hierarchical_info), r_lvl, ctx$hierarchical_info, allow_identity = TRUE)) {
                local_struct_count <- local_struct_count + 1
              }
            }
            
            use_partitioning <- (local_struct_count == 1) && identical(s_lvl, r_lvl)

            s_def <- jags_structure_definition(
              ctx$structures[[s_name]],
              variable_name = response,
              s_name = s_name,
              loop_bound = s_bound,
              zeros_name = s_zeros,
              is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
              i_index = s_idx_var,
              use_partitioning = use_partitioning,
              engine = ctx$engine
            )

            ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
            
            if (isTRUE(s_def$partition_handled)) {
              prec_nm     <- paste0("tau_u_", s_name, "_", response)
              tau_res_nm  <- paste0("tau_res_", response)
              sigma_ph_nm <- paste0("sigma_", s_name, "_", response)
              sigma_rs_nm <- paste0("sigma_", response, "_res")
              sigma_tot_nm<- paste0("sigma_total_", response)
              
              nodes_to_add <- c(prec_nm, tau_res_nm, sigma_ph_nm, sigma_rs_nm, sigma_tot_nm)
              if (is.environment(ctx$declared_nodes)) {
                ctx$declared_nodes$nodes <- unique(c(ctx$declared_nodes$nodes, nodes_to_add))
              }
              ctx$partitioned_responses <- unique(c(ctx$partitioned_responses, response))
            }

            if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
               additive_terms <- paste0(additive_terms, " + ", s_def$term)
            }

            param_name <- paste0("sigma_", response, "_", s_name)
            tau_name <- paste0("tau_u_", response, "_", s_name)
            
            existing_idx <- which(sapply(ctx$param_map, function(p) {
              p$response == response && p$predictor == s_name && p$type == "structure"
            }))
            
            is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(u, s_name), logical(1)))
            
            if (is_unified) {
              unified_tau <- paste0("tau_u_", s_name, "_", response)
              unified_sigma <- paste0("sigma_", s_name, "_", response)
              ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
              ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
              
              if (isTRUE(s_def$partition_handled)) {
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(
                  response = response, 
                  predictor = s_name, 
                  parameter = paste0("lambda_", response), 
                  type = "structure"
                )
              }
            } else {
              if (length(existing_idx) > 0) {
                existing_param <- ctx$param_map[[existing_idx[1]]]$parameter
                if (grepl(paste0("_", s_name, "$"), param_name) && !grepl(paste0("_", s_name, "$"), existing_param)) {
                   ctx$param_map[[existing_idx[1]]]$parameter <- tau_name
                   ctx$param_map[[existing_idx[2]]]$parameter <- param_name
                }
              } else {
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = tau_name, equation_index = NA, type = "structure")
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = param_name, equation_index = NA, type = "structure")
              }
            }
          }

          # 2. Random Effects (Grouped)
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
            
            if (r_name %in% names(ctx$structures) || any(vapply(names(ctx$structures), function(s) grepl(s, r_name), logical(1)))) {
              is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
              if (is_unified_r) next
            }
            
            s_suffix_r <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix_r)
            u <- paste0("u_", response, suffix, s_suffix_r)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix_r)
            
            r_lvl <- get_random_level(response, r_name, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
            n_groups <- paste0("N_", r_name)
            zeros_name <- paste0("zeros_", r_name)
            prec_name <- paste0("Prec_", r_name)

            ctx$model_lines <- safe_add_lines(ctx$model_lines, c(
              paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
              paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- (", u_std, "[g] - mean(", u_std, "[1:", n_groups, "])) / sqrt(", tau_u, ") }")
            ), ctx$declared_nodes)
            
            group_idx <- get_group_idx_string(response, r_name, ctx$hierarchical_info, default_N = ctx$main_loop_N)
            additive_terms <- paste0(additive_terms, " + ", u, "[", group_idx, "]")
            
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = r_name, parameter = tau_u, equation_index = NA, type = "structure")
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name), equation_index = NA, type = "structure")
          }

          # 3. Unified Observation Loop
          tau_res <- paste0("tau_res_", response, suffix)
          
          partition_already_handled <- FALSE
          if (!is.null(ctx$structures)) {
            if (exists("s_def") && isTRUE(s_def$partition_handled) && s_def$variable_name == response) {
              partition_already_handled <- TRUE
            }
          }
          
          if (!partition_already_handled) {
            p_obj <- get_family_object(ctx$dist_list[[response]] %||% "gaussian")
            ctx$model_lines <- safe_add_lines(ctx$model_lines, jags_family_precision_prior(p_obj, tau_res), ctx$declared_nodes)
          }

          ctx$model_lines <- c(
            ctx$model_lines,
            paste0("  for (i in 1:", loop_bound, ") {"),
            paste0(
              "    ",
              response_var,
              "[i] ~ dnorm(",
              mu,
              "[i]",
              additive_terms,
              ", ",
              tau_res,
              ")"
            )
          )

          if (ctx$engine == "jags") {
            ctx$model_lines <- c(
              ctx$model_lines,
              paste0(
                "    log_lik_",
                response,
                suffix,
                "[i] <- logdensity.norm(",
                response_var,
                "[i], ",
                mu,
                "[i]",
                additive_terms,
                ", ",
                tau_res,
                ")"
              )
            )
          }
          ctx$model_lines <- safe_add_lines(ctx$model_lines, "  }", ctx$declared_nodes)
        }
      } else if (dist == "binomial") {
        err <- paste0("err_", response, suffix)

        if (ctx$independent) {
          ctx$model_lines <- c(
            ctx$model_lines,
            paste0("  # Independent (Standard GLM) for binomial: ", response)
          )
        } else {
          epsilon <- paste0("epsilon_", response, suffix)
          tau_res <- paste0("tau_res_", response, suffix)
          
          total_u <- ""
          local_obs_code <- c()

          for (s_name in ctx$structure_names) {
            if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next
            
            s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
            s_bound <- if (is.null(s_lvl)) loop_bound else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

            s_def <- jags_structure_definition(
              ctx$structures[[s_name]],
              variable_name = response,
              s_name = s_name,
              loop_bound = s_bound,
              zeros_name = s_zeros,
              category_index = NULL,
              is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
              i_index = s_idx_var,
              engine = ctx$engine
            )
            is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(u, s_name), logical(1)))
            if (!is.null(s_def)) {
              if (is_unified) {
                unified_tau <- paste0("tau_u_", s_name, "_", response)
                unified_sigma <- paste0("sigma_", s_name, "_", response)
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
              } else {
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
              }

              ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
              
              if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
                total_u <- paste0(total_u, " + ", s_def$term)
              }
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
            
            if (r_name %in% names(ctx$structures) || any(vapply(names(ctx$structures), function(s) grepl(s, r_name), logical(1)))) {
              is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
              if (is_unified_r) next
            }

            s_suffix_r <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix_r)
            u <- paste0("u_", response, suffix, s_suffix_r)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix_r)
            
            r_lvl <- get_random_level(response, r_name, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
            n_groups <- paste0("N_", r_name)
            zeros_name <- paste0("zeros_", r_name)
            prec_name <- paste0("Prec_", r_name)

            ctx$model_lines <- safe_add_lines(ctx$model_lines, c(
              paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
              paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- (", u_std, "[g] - mean(", u_std, "[1:", n_groups, "])) / sqrt(", tau_u, ") }")
            ), ctx$declared_nodes)
            
            group_idx <- get_group_idx_string(response, r_name, ctx$hierarchical_info, default_N = ctx$main_loop_N)
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
            
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
      } else if (dist == "multinomial") {
        K_var <- paste0("K_", response)
        err <- paste0("err_", response, suffix)

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
              s_obj <- ctx$structures[[s_idx]]

              if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next

              s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
              s_bound <- if (is.null(s_lvl)) get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N) else paste0("N_", s_lvl)
              s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
              s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

              s_def <- jags_structure_definition(
                s_obj,
                variable_name = response,
                s_name = s_name,
                loop_bound = s_bound,
                zeros_name = s_zeros,
                category_index = "k",
                is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
                i_index = s_idx_var,
                engine = ctx$engine
              )

              if (!is.null(s_def)) {
                ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
                total_u <- paste0(total_u, " + ", s_def$term)
                
                is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
                if (is_unified) {
                  unified_tau <- paste0("tau_u_", s_name, "_", response, "[k]")
                  unified_sigma <- paste0("sigma_", s_name, "_", response, "[k]")
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
                } else {
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name, "[k]"), equation_index = NA, type = "structure")
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name, "[k]"), equation_index = NA, type = "structure")
                }
              }
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
            
            if (r_name %in% names(ctx$structures) || any(vapply(names(ctx$structures), function(s) grepl(s, r_name), logical(1)))) {
              is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
              if (is_unified_r) next
            }

            s_suffix_r <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix_r)
            u <- paste0("u_", response, suffix, s_suffix_r)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix_r)
            
            r_lvl <- get_random_level(response, r_name, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
            n_groups <- paste0("N_", r_name)
            zeros_name <- paste0("zeros_", r_name)
            prec_name <- paste0("Prec_", r_name)

            ctx$model_lines <- safe_add_lines(ctx$model_lines, c(
              paste0("  ", u_std, "[1:", n_groups, ", k] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
              paste0("  for (g in 1:", n_groups, ") { ", u, "[g, k] <- (", u_std, "[g, k] - mean(", u_std, "[1:", n_groups, ", k])) / sqrt(", tau_u, "[k]) }")
            ), ctx$declared_nodes)
            
            group_idx <- get_group_idx_string(response, r_name, ctx$hierarchical_info, default_N = ctx$main_loop_N)
            total_u <- paste0(total_u, " + ", u, "[", group_idx, ", k]")
            
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0(tau_u, "[k]"), equation_index = NA, type = "structure")
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
      } else if (dist == "ordinal") {
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
              s_obj <- ctx$structures[[s_idx]]

              if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next

              s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
              s_bound <- if (is.null(s_lvl)) get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N) else paste0("N_", s_lvl)
              s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
              s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

              s_def <- jags_structure_definition(
                s_obj,
                variable_name = response,
                s_name = s_name,
                loop_bound = s_bound,
                zeros_name = s_zeros,
                is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
                i_index = s_idx_var,
                engine = ctx$engine
              )

              if (!is.null(s_def)) {
                ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
                total_u <- paste0(total_u, " + ", s_def$term)

                is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
                if (is_unified) {
                  unified_tau <- paste0("tau_u_", s_name, "_", response)
                  unified_sigma <- paste0("sigma_", s_name, "_", response)
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
                } else {
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
                }
              }
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
            
            if (r_name %in% names(ctx$structures) || any(vapply(names(ctx$structures), function(s) grepl(s, r_name), logical(1)))) {
              is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
              if (is_unified_r) next
            }

            s_suffix_r <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix_r)
            u <- paste0("u_", response, suffix, s_suffix_r)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix_r)
            
            r_lvl <- get_random_level(response, r_name, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
            n_groups <- paste0("N_", r_name)
            zeros_name <- paste0("zeros_", r_name)
            prec_name <- paste0("Prec_", r_name)

            ctx$model_lines <- safe_add_lines(ctx$model_lines, c(
              paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
              paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- (", u_std, "[g] - mean(", u_std, "[1:", n_groups, "])) / sqrt(", tau_u, ") }")
            ), ctx$declared_nodes)
            
            group_idx <- get_group_idx_string(response, r_name, ctx$hierarchical_info, default_N = ctx$main_loop_N)
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
            
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
      } else if (dist == "poisson") {
        err <- paste0("err_", response, suffix)

        if (ctx$independent) {
          ctx$model_lines <- c(
            ctx$model_lines,
            paste0("  # Independent Poisson (no random effect term)"),
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
              s_obj <- ctx$structures[[s_idx]]

              if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next

              s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
              s_bound <- if (is.null(s_lvl)) get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N) else paste0("N_", s_lvl)
              s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
              s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

              s_def <- jags_structure_definition(
                s_obj,
                variable_name = response,
                s_name = s_name,
                loop_bound = s_bound,
                zeros_name = s_zeros,
                is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
                i_index = s_idx_var,
                engine = ctx$engine
              )

              if (!is.null(s_def)) {
                ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
                total_u <- paste0(total_u, " + ", s_def$term)

                is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
                if (is_unified) {
                  unified_tau <- paste0("tau_u_", s_name, "_", response)
                  unified_sigma <- paste0("sigma_", s_name, "_", response)
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
                } else {
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                  ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
                }
              }
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
            
            if (r_name %in% names(ctx$structures) || any(vapply(names(ctx$structures), function(s) grepl(s, r_name), logical(1)))) {
              is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
              if (is_unified_r) next
            }

            s_suffix_r <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix_r)
            u <- paste0("u_", response, suffix, s_suffix_r)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix_r)
            
            r_lvl <- get_random_level(response, r_name, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
            n_groups <- paste0("N_", r_name)
            zeros_name <- paste0("zeros_", r_name)
            prec_name <- paste0("Prec_", r_name)

            ctx$model_lines <- safe_add_lines(ctx$model_lines, c(
              paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
              paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- (", u_std, "[g] - mean(", u_std, "[1:", n_groups, "])) / sqrt(", tau_u, ") }")
            ), ctx$declared_nodes)
            
            group_idx <- get_group_idx_string(response, r_name, ctx$hierarchical_info, default_N = ctx$main_loop_N)
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
            
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
      } else if (dist == "occupancy") {
        total_u <- ""
        if (length(ctx$structures) > 0) {
          for (s_idx in seq_along(ctx$structures)) {
            s_name <- names(ctx$structures)[s_idx]
            if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
            s_obj <- ctx$structures[[s_idx]]

            if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next

            s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
            s_bound <- if (is.null(s_lvl)) get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N) else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

            s_def <- jags_structure_definition(
              s_obj,
              variable_name = response,
              s_name = s_name,
              loop_bound = s_bound,
              zeros_name = s_zeros,
              is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
              i_index = s_idx_var,
              engine = ctx$engine
            )

            if (!is.null(s_def)) {
              ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
              total_u <- paste0(total_u, " + ", s_def$term)

              is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
              if (is_unified) {
                unified_tau <- paste0("tau_u_", s_name, "_", response)
                unified_sigma <- paste0("sigma_", s_name, "_", response)
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
              } else {
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
              }
            }
          }
        }

        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("  # Occupancy Model for ", response),
          paste0("  for (i in 1:", eq_loop_N, ") {"),
          paste0("    logit(psi_", response, "[i]) <- mu_", response, suffix, "[i]", total_u)
        )

        fam_obj <- get_family_object(dist)
        curr_pred <- NULL
        for (eq in ctx$eq_list) {
          if (eq$response == response) {
            curr_pred <- eq$predictors
            break
          }
        }

        def <- jags_family_definition(fam_obj, response, curr_pred)
        if (!is.null(def$model_code)) {
          ctx$model_lines <- safe_add_lines(ctx$model_lines, def$model_code, ctx$declared_nodes)
        } else {
          stop(paste("Unknown distribution or missing module for:", dist))
        }

        ctx$model_lines <- c(ctx$model_lines, "  }")
      } else if (dist == "negbinomial" || dist == "zinb") {
        err <- paste0("err_", response, suffix)

        epsilon <- paste0("epsilon_", response, suffix)
        tau_res <- paste0("tau_res_", response, suffix)

        total_u <- ""

        if (length(ctx$structures) > 0) {
          for (s_idx in seq_along(ctx$structures)) {
            s_name <- names(ctx$structures)[s_idx]
            if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
            s_obj <- ctx$structures[[s_idx]]

            if (!is_valid_structure_mapping(get_struct_lvl(s_name, ctx$hierarchical_info), get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) next

            s_lvl <- get_struct_lvl(s_name, ctx$hierarchical_info)
            s_bound <- if (is.null(s_lvl)) get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N) else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_idx_var <- get_struct_index(s_name, response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

            s_def <- jags_structure_definition(
              s_obj,
              variable_name = response,
              s_name = s_name,
              loop_bound = s_bound,
              zeros_name = s_zeros,
              is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
              i_index = s_idx_var,
              engine = ctx$engine
            )

            if (!is.null(s_def)) {
              ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
              total_u <- paste0(total_u, " + ", s_def$term)

              is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
              if (is_unified) {
                unified_tau <- paste0("tau_u_", s_name, "_", response)
                unified_sigma <- paste0("sigma_", s_name, "_", response)
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
              } else {
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
              }
            }
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
          
          if (r_name %in% names(ctx$structures) || any(vapply(names(ctx$structures), function(s) grepl(s, r_name), logical(1)))) {
            is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
            if (is_unified_r) next
          }

          s_suffix_r <- paste0("_", r_name)
          u_std <- paste0("u_std_", response, suffix, s_suffix_r)
          u <- paste0("u_", response, suffix, s_suffix_r)
          tau_u <- paste0("tau_u_", response, suffix, s_suffix_r)
          
          r_lvl <- get_random_level(response, r_name, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
          n_groups <- paste0("N_", r_name)
          zeros_name <- paste0("zeros_", r_name)
          prec_name <- paste0("Prec_", r_name)

          ctx$model_lines <- safe_add_lines(ctx$model_lines, c(
            paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
            paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- (", u_std, "[g] - mean(", u_std, "[1:", n_groups, "])) / sqrt(", tau_u, ") }")
          ), ctx$declared_nodes)
          
          group_idx <- get_group_idx_string(response, r_name, ctx$hierarchical_info, default_N = ctx$main_loop_N)
          total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
          
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
    }
  }

  # Likelihoods for MAG correlated variables
  for (var in names(ctx$vars_error_terms)) {
    dist <- if (!is.null(ctx$dist_list[[var]])) ctx$dist_list[[var]] else "gaussian"
    if (dist != "gaussian") next

    err_terms <- ctx$vars_error_terms[[var]]
    suffix <- if ((ctx$response_counter[[var]] %||% 0) > 1) "1" else ""
    loop_bound <- get_loop_bound(var, ctx$hierarchical_info, default_N = ctx$main_loop_N)

    structure_term_str <- ""
    structure_terms <- character()
    if (length(ctx$structures) > 0) {
      for (s_idx in seq_along(ctx$structures)) {
        s_name <- names(ctx$structures)[s_idx]
        if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
        s_obj <- ctx$structures[[s_idx]]

        local_struct_count <- 0
        for (sn in names(ctx$structures)) {
          if (is_valid_structure_mapping(get_struct_lvl(sn, ctx$hierarchical_info), get_var_level(var, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars), ctx$hierarchical_info, allow_identity = TRUE)) {
            local_struct_count <- local_struct_count + 1
          }
        }
        
        use_partitioning_flag <- (dist == "gaussian" && local_struct_count == 1 && length(err_terms) == 0)

        s_def <- jags_structure_definition(
          s_obj,
          variable_name = var,
          s_name = s_name,
          loop_bound = loop_bound,
          is_multi = is_struct_multi(s_name, ctx$hierarchical_info),
          use_partitioning = use_partitioning_flag,
          engine = ctx$engine
        )

        if (!is.null(s_def$model_lines)) {
          ctx$model_lines <- safe_add_lines(ctx$model_lines, s_def$model_lines, ctx$declared_nodes)
        }
        if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
          structure_terms <- c(structure_terms, s_def$term)
        }

        tau_u_var <- paste0("tau_u_", var, "_", s_name)
        sigma_u_var <- paste0("sigma_", var, "_", s_name)
        ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = var, predictor = s_name, parameter = tau_u_var, equation_index = NA, type = "structure")
        ctx$param_map[[length(ctx$param_map) + 1]] <- list(response = var, predictor = s_name, parameter = sigma_u_var, equation_index = NA, type = "structure")
      }
      if (length(structure_terms) > 0) {
        structure_term_str <- paste0(" + ", paste(structure_terms, collapse = " + "))
      }
    }

    error_sum <- paste(c(err_terms), collapse = " + ")

    mu_var <- paste0("mu_", var, suffix)
    tau_res_var <- paste0("tau_res_", var, suffix)

    ctx$model_lines <- c(
      ctx$model_lines,
      paste0("  # Likelihood for ", var, " (with correlated residual errors)"),
      paste0("  for (i in 1:", loop_bound, ") {"),
      paste0(
        "    ",
        var,
        "[i] ~ dnorm(",
        mu_var,
        "[i] + ",
        error_sum,
        structure_term_str,
        ", ",
        tau_res_var,
        ")"
      )
    )
    if (ctx$engine == "jags") {
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0(
          "    log_lik_",
          var,
          suffix,
          "[i] <- logdensity.norm(",
          var,
          "[i], ",
          mu_var,
          "[i] + ",
          error_sum,
          structure_term_str,
          ", ",
          tau_res_var,
          ")"
        )
      )
    }
    ctx$model_lines <- safe_add_lines(ctx$model_lines, "  }", ctx$declared_nodes)
  }

  invisible(ctx)
}
