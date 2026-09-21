#' Generate linear predictors for JAGS/BUGS model
#'
#' Constructs the systematic linear predictors (mu[i] <- ...) across equations
#' for all supported distributions (Gaussian, Binomial, Poisson, Multinomial,
#' Ordinal, Negative Binomial, ZIP, ZINB, Occupancy).
#'
#' @param ctx Builder context environment
#' @keywords internal
generate_linear_predictors <- function(ctx) {
  for (j in seq_along(ctx$eq_list)) {
    eq <- ctx$eq_list[[j]]
    response <- eq$response
    predictors <- eq$predictors
    dist <- ctx$dist_list[[response]] %||% "gaussian"

    # Hierarchical Loop Wrapper
    eq_loop_N <- get_loop_bound(response, ctx$hierarchical_info, default_N = ctx$main_loop_N)
    ctx$model_lines <- safe_add_lines(ctx$model_lines, paste0("  for (i in 1:", eq_loop_N, ") {"), ctx$declared_nodes)

    # Check if this is a deterministic identity equation
    is_identity <- FALSE
    if (length(predictors) == 1) {
      if (!is.null(ctx$categorical_vars)) {
        is_dummy <- any(sapply(ctx$categorical_vars, function(cv) {
          response %in% cv$dummies
        }))
        if (is_dummy) is_identity <- TRUE
      }
    }

    if (is_identity) {
      expr <- term_to_jags_expression(predictors[1])
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    ", response, "[i] <- ", expr),
        "  }"
      )
      next
    }

    # Count and assign unique suffix for the response variable
    response_count <- ctx$response_counter[[response]] %||% 0
    response_count <- response_count + 1
    ctx$response_counter[[response]] <- response_count
    suffix <- if (response_count == 1) "" else as.character(response_count)

    alpha <- paste0("alpha_", response, suffix)
    ctx$param_map[[length(ctx$param_map) + 1]] <- list(
      response = response,
      predictor = "(Intercept)",
      parameter = alpha,
      equation_index = j,
      type = "coefficient"
    )
    linpred <- alpha
    for (pred in predictors) {
      key <- paste(response, pred, suffix, sep = "_")
      if (!key %in% names(ctx$beta_counter)) {
        beta_name <- paste0("beta_", response, suffix, "_", pred)
        ctx$beta_counter[[key]] <- beta_name
      }
      beta_name <- ctx$beta_counter[[key]]

      pred_dist <- ctx$dist_list[[pred]] %||% "gaussian"
      resp_level <- get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
      pred_idx <- get_pred_index(pred, resp_level, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

      pred_fam_obj <- get_family_object(pred_dist)
      actual_pred_idx <- get_latent_predictor_name(pred_fam_obj, pred, pred_idx)
      linpred <- paste0(linpred, " + ", beta_name, "*", actual_pred_idx)

      ctx$param_map[[length(ctx$param_map) + 1]] <- list(
        response = response,
        predictor = pred,
        parameter = beta_name,
        equation_index = j,
        type = "coefficient"
      )
    }

    if (dist == "gaussian" || dist == "occupancy" || grepl("^p_", response)) {
      mu <- paste0("mu_", response, suffix)
      ctx$model_lines <- safe_add_lines(ctx$model_lines, paste0("    ", mu, "[i] <- ", linpred), ctx$declared_nodes)
    } else if (dist == "binomial") {
      mu_err <- paste0("mu_err_", response, suffix)
      err_name <- paste0("err_", response, suffix)
      p <- paste0("p_", response, suffix)
      
      err_term <- ""
      if (!ctx$independent || err_name %in% names(ctx$vars_error_terms)) {
        err_term <- paste0(" + ", err_name, "[i]")
      }

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    ", mu_err, "[i] <- 0"),
        paste0("    logit(", p, "[i]) <- ", linpred, err_term),
        paste0("    ", response, "[i] ~ dbern(", p, "[i])")
      )
      if (ctx$engine == "jags") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    log_lik_", response, suffix, "[i] <- logdensity.bern(", response, "[i], ", p, "[i])")
        )
      }
    } else if (dist == "multinomial") {
      K_var <- paste0("K_", response)
      err <- paste0("err_", response, suffix)

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    # Multinomial linear predictor for ", response),
        paste0("    L_", response, "[i, 1] <- 0"),
        paste0("    for (k in 2:", K_var, ") {")
      )

      linpred_k <- paste0("alpha_", response, "[k]")

      for (pred in predictors) {
        beta_name <- paste0("beta_", response, "_", pred)
        resp_level <- get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
        pred_idx <- get_pred_index(pred, resp_level, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
        linpred_k <- paste0(linpred_k, " + ", beta_name, "[k] * ", pred_idx)

        if (!paste0("alpha_", response) %in% names(ctx$beta_counter)) {
          ctx$beta_counter[[paste0("alpha_", response)]] <- TRUE
          ctx$param_map[[length(ctx$param_map) + 1]] <- list(
            response = response,
            predictor = "(Intercepts)",
            parameter = paste0("alpha_", response, "[]"),
            equation_index = j,
            type = "coefficient"
          )
        }

        key <- paste(response, pred, suffix, sep = "_")
        if (!key %in% names(ctx$beta_counter)) {
          ctx$beta_counter[[key]] <- beta_name
          ctx$param_map[[length(ctx$param_map) + 1]] <- list(
            response = response,
            predictor = pred,
            parameter = paste0(beta_name, "[]"),
            equation_index = j,
            type = "coefficient"
          )
        }
      }

      if (is.null(linpred_k) || nchar(trimws(linpred_k)) == 0) linpred_k <- "0"
      linpred_k <- paste0(linpred_k, " + ", err, "[i, k]")

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("      L_", response, "[i, k] <- max(-20, min(20, ", linpred_k, "))"),
        "    }",
        paste0("    # Softmax for ", response)
      )

      if (ctx$engine == "nimble") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    max_L_", response, "[i] <- max(L_", response, "[i, 1:", K_var, "])"),
          paste0("    for (k in 1:", K_var, ") {"),
          paste0("      exp_L_", response, "[i, k] <- exp(L_", response, "[i, k] - max_L_", response, "[i])"),
          "    }"
        )
      } else {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    for (k in 1:", K_var, ") {"),
          paste0("      exp_L_", response, "[i, k] <- exp(L_", response, "[i, k])"),
          "    }"
        )
      }

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    sum_exp_L_", response, "[i] <- sum(exp_L_", response, "[i, 1:", K_var, "])"),
        paste0("    for (k in 1:", K_var, ") {"),
        paste0("      p_", response, "[i, k] <- (exp_L_", response, "[i, k] / sum_exp_L_", response, "[i]) * 0.9999999 + 1.0e-10"),
        "    }",
        paste0("    ", response, "[i] ~ dcat(p_", response, "[i, 1:", K_var, "])")
      )
      if (ctx$engine == "jags") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    log_lik_", response, suffix, "[i] <- logdensity.cat(", response, "[i], p_", response, "[i, 1:", K_var, "])")
        )
      }
    } else if (dist == "ordinal") {
      K_var <- paste0("K_", response)
      err_name <- paste0("err_", response, suffix)
      eta <- paste0("eta_", response, suffix)

      err_term <- ""
      if (!ctx$independent || err_name %in% names(ctx$vars_error_terms)) {
        err_term <- paste0(" + ", err_name, "[i]")
      }

      linpred_no_int <- "0"
      if (!is.null(predictors) && length(predictors) > 0) {
        for (pred in predictors) {
          beta_name <- paste0("beta_", response, "_", pred)
          resp_level <- get_var_level(response, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)
          pred_idx <- get_pred_index(pred, resp_level, ctx$hierarchical_info, equations = ctx$equations, latent = ctx$latent, categorical_vars = ctx$categorical_vars)

          linpred_no_int <- paste0(linpred_no_int, " + ", beta_name, " * ", pred_idx)

          key <- paste(response, pred, suffix, sep = "_")
          if (!beta_name %in% names(ctx$beta_counter)) {
            ctx$beta_counter[[length(ctx$beta_counter) + 1]] <- beta_name
            names(ctx$beta_counter)[length(ctx$beta_counter)] <- key
            ctx$param_map[[length(ctx$param_map) + 1]] <- list(
              response = response,
              predictor = pred,
              parameter = beta_name,
              equation_index = j,
              type = "coefficient"
            )
          }
        }
      }

      ctx$param_map[[length(ctx$param_map) + 1]] <- list(
        response = response,
        predictor = "(Cutpoints)",
        parameter = paste0("cutpoint_", response, suffix),
        equation_index = j,
        type = "coefficient"
      )

      if (is.null(err_term) || nchar(trimws(err_term)) == 0) err_term <- ""

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    # Ordinal linear predictor for ", response),
        paste0("    ", eta, "[i] <- ", linpred_no_int, err_term),
        paste0("    for (k in 1:(", K_var, "-1)) {"),
        paste0("      logit(Q_", response, "[i, k]) <- cutpoint_", response, "[k] - ", eta, "[i]"),
        "    }",
        paste0("    p_", response, "[i, 1] <- Q_", response, "[i, 1]"),
        paste0("    for (k in 2:(", K_var, "-1)) {"),
        paste0("      p_", response, "[i, k] <- Q_", response, "[i, k] - Q_", response, "[i, k-1]"),
        "    }",
        paste0("    p_", response, "[i, ", K_var, "] <- 1 - Q_", response, "[i, ", K_var, "-1]"),
        paste0("    ", response, "[i] ~ dcat(p_", response, "[i, 1:", K_var, "])")
      )
      if (ctx$engine == "jags") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    log_lik_", response, suffix, "[i] <- logdensity.cat(", response, "[i], p_", response, "[i, 1:", K_var, "])")
        )
      }
    } else if (dist == "poisson") {
      mu <- paste0("mu_", response, suffix)
      err_term <- ""
      err_name <- paste0("err_", response, suffix)
      
      if (!ctx$independent || err_name %in% names(ctx$vars_error_terms)) {
        err_term <- paste0(" + ", err_name, "[i]")
      }

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    # Poisson log link for ", response)
      )
      
      if (ctx$engine == "jags") {
        ctx$model_lines <- c(ctx$model_lines, paste0(
          "    ", mu, "[i] <- exp(max(-20, min(10, ", linpred, err_term, ")))"
        ))
      } else {
        ctx$model_lines <- c(ctx$model_lines, paste0(
          "    ", mu, "[i] <- exp(", linpred, err_term, ")"
        ))
      }
      
      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    ", response, "[i] ~ dpois(", mu, "[i])")
      )
      if (ctx$engine == "jags") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    log_lik_", response, suffix, "[i] <- logdensity.pois(", response, "[i], ", mu, "[i])")
        )
      }
    } else if (dist == "negbinomial") {
      mu <- paste0("mu_", response, suffix)
      p <- paste0("p_", response, suffix)
      r <- paste0("r_", response, suffix)
      
      err_term <- ""
      err_name <- paste0("err_", response, suffix)
      
      if (!ctx$independent || err_name %in% names(ctx$vars_error_terms)) {
        err_term <- paste0(" + ", err_name, "[i]")
      }

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    # Negative Binomial log link for ", response)
      )
      
      if (ctx$engine == "jags") {
        ctx$model_lines <- c(ctx$model_lines, paste0(
          "    ", mu, "[i] <- exp(max(-20, min(10, ", linpred, err_term, ")))"
        ))
      } else {
        ctx$model_lines <- c(ctx$model_lines, paste0(
          "    ", mu, "[i] <- exp(", linpred, err_term, ")"
        ))
      }
      if (ctx$engine == "jags") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    ", p, "[i] <- ", r, " / (", r, " + ", mu, "[i])"),
          paste0("    ", response, "[i] ~ dnegbin(", p, "[i], ", r, ")")
        )
      } else {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    ", response, "[i] ~ dnb_because(", mu, "[i], ", r, ")")
        )
      }
      if (ctx$engine == "jags") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    log_lik_", response, suffix, "[i] <- logdensity.negbin(", response, "[i], ", p, "[i], ", r, ")")
        )
      }
    } else if (dist == "zip") {
      err <- paste0("err_", response, suffix)
      mu <- paste0("mu_", response, suffix)
      psi <- paste0("psi_", response, suffix)

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    # ZIP log link for ", response),
        paste0("    ", mu, "[i] <- exp(max(-20, min(10, ", linpred, " + ", err, "[i])))"),
        paste0("    lik_zero_", response, "[i] <- ", psi, " + (1-", psi, ") * exp(-", mu, "[i])"),
        paste0(
          "    lik_pos_", response, "[i] <- (1-", psi, ") * exp(",
          if (ctx$engine == "jags") {
            paste0("logdensity.pois(", response, "[i], ", mu, "[i])")
          } else {
            paste0("dpois(", response, "[i], ", mu, "[i], 1)")
          },
          ")"
        )
      )

      if (ctx$engine == "jags") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    is_zero_", response, suffix, "[i] <- step(0.5 - ", response, "[i])"),
          paste0("    lik_", response, suffix, "[i] <- is_zero_", response, suffix, "[i] * lik_zero_", response, suffix, "[i] + (1 - is_zero_", response, suffix, "[i]) * lik_pos_", response, suffix, "[i]"),
          paste0("    log_lik_", response, suffix, "[i] <- log(max(1.0E-30, lik_", response, suffix, "[i]))"),
          paste0("    phi_", response, suffix, "[i] <- -log_lik_", response, suffix, "[i] + 10000"),
          paste0("    zeros[i] ~ dpois(phi_", response, suffix, "[i])")
        )
      } else {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    ", response, "[i] ~ dzip_because(", mu, "[i], psi_", response, suffix, ")")
        )
      }
    } else if (dist == "zinb") {
      err <- paste0("err_", response, suffix)
      mu <- paste0("mu_", response, suffix)
      r <- paste0("r_", response, suffix)
      p <- paste0("p_", response, suffix)
      psi <- paste0("psi_", response, suffix)

      ctx$model_lines <- c(
        ctx$model_lines,
        paste0("    # ZINB log link for ", response),
        paste0("    ", mu, "[i] <- exp(max(-20, min(10, ", linpred, " + ", err, "[i])))"),
        paste0("    ", p, "[i] <- ", r, " / (", r, " + ", mu, "[i])")
      )

      if (ctx$engine == "jags") {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    lik_zero_", response, suffix, "[i] <- ", psi, " + (1-", psi, ") * pow(", p, "[i], ", r, ")"),
          paste0(
            "    lik_pos_", response, suffix, "[i] <- (1-", psi, ") * exp(",
            if (ctx$engine == "jags") {
              paste0("logdensity.negbin(", response, "[i], ", p, "[i], ", r, ")")
            } else {
              paste0("dnegbin(", response, "[i], ", p, "[i], ", r, ", 1)")
            },
            ")"
          ),
          paste0("    is_zero_", response, suffix, "[i] <- step(0.5 - ", response, "[i])"),
          paste0("    lik_", response, suffix, "[i] <- is_zero_", response, suffix, "[i] * lik_zero_", response, suffix, "[i] + (1 - is_zero_", response, suffix, "[i]) * lik_pos_", response, suffix, "[i]"),
          paste0("    log_lik_zero_", response, suffix, "[i] <- log(max(1.0E-30, lik_zero_", response, suffix, "[i]))"),
          paste0(
            "    log_lik_pos_", response, suffix, "[i] <- log(max(1.0E-30, 1-", psi, ")) + ",
            if (ctx$engine == "jags") {
              paste0("logdensity.negbin(", response, "[i], ", p, "[i], ", r, ")")
            } else {
              paste0("dnegbin(", response, "[i], ", p, "[i], ", r, ", 1)")
            }
          ),
          paste0("    log_lik_", response, suffix, "[i] <- is_zero_", response, suffix, "[i] * log_lik_zero_", response, suffix, "[i] + (1 - is_zero_", response, suffix, "[i]) * log_lik_pos_", response, suffix, "[i]"),
          paste0("    phi_", response, suffix, "[i] <- -log_lik_", response, suffix, "[i] + 10000"),
          paste0("    zeros_", response, "[i] ~ dpois(phi_", response, suffix, "[i])")
        )
      } else {
        ctx$model_lines <- c(
          ctx$model_lines,
          paste0("    ", response, "[i] ~ dzinb_because(", mu, "[i], ", r, ", psi_", response, suffix, ")")
        )
      }
    } else {
      stop(paste("Unknown distribution:", dist))
    }

    ctx$model_lines <- c(ctx$model_lines, "  }")
  }

  invisible(ctx)
}
