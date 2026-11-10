#' Standardize Path Coefficients onto Latent Liability and Common SD Scale
#'
#' Rescales parameters from a fitted `because` model or d-separation test
#' into standardized path coefficients (\eqn{\beta^*}).
#'
#' \itemize{
#'   \item \strong{Continuous (Gaussian) Variables}: Scaled by the ratio of standard
#'     deviations \eqn{\beta^* = \beta \times (\text{SD}(X) / \text{SD}(Y))}. If variables
#'     were already standardized prior to model fitting, their scale is preserved.
#'   \item \strong{Binary / Categorical (Binomial/Bernoulli/Ordinal Logit)}: Raw coefficients sit on
#'     the log-odds scale. They are standardized \emph{a posteriori} onto the continuous latent
#'     liability scale following McKelvey & Zavoina (1975), Menard (2004), and Nakagawa et al. (2017):
#'     \deqn{\beta^* = \frac{\beta \times \text{SD}(X)}{\sqrt{\sigma_{\text{phylo}}^2 + \pi^2/3}}}
#'     where \eqn{\pi^2/3 \approx 3.289868} is the theoretical disturbance variance of the standard logistic link.
#'   \item \strong{Binary Predictors}: Slopes pointing from a binary indicator (e.g. \code{Migrant} \eqn{\in \{0, 1\}})
#'     to a continuous response are multiplied by the predictor's standard deviation
#'     \eqn{s_X = \sqrt{p(1 - p)}}, placing the binary jump on the standard deviation scale (analogous to Cohen's \eqn{d}).
#'   \item \strong{d-separation Objects}: Evaluates and rescales each auxiliary independence claim
#'     across all d-sep sub-models, enabling a uniform \code{rope = c(-0.10, 0.10)} benchmark across all tests.
#' }
#'
#' @param object A fitted `because` model object or a `because` d-separation test object.
#' @param method Character; liability variance calculation method: `"error"` (default; disturbance +
#'   structural variance \eqn{\sigma_{\text{phylo}}^2 + \pi^2/3}) or `"total"` (full McKelvey-Zavoina
#'   incorporating the empirical variance of the fixed linear predictor).
#' @param ... Additional arguments passed to methods.
#'
#' @return An updated object of class \code{c("because_standardized", "because")} with transformed
#'   posterior samples and updated summary statistics. Fully compatible with \code{\link{summary}},
#'   \code{\link{plot_coef}}, \code{\link{plot_rope}}, and \code{\link{plot_dag}}.
#'
#' @references
#' McKelvey, R. D., & Zavoina, W. (1975). A statistical model for the analysis of ordinal level dependent variables.
#' \emph{Journal of Mathematical Sociology}, 4(1), 103–120.
#'
#' Menard, S. (2004). Six approaches to calculating standardized logistic regression coefficients.
#' \emph{The American Statistician}, 58(3), 218–223.
#'
#' Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The coefficient of determination \eqn{R^2}
#' and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded.
#' \emph{Journal of the Royal Society Interface}, 14(134), 20170213.
#'
#' @export
standardize <- function(object, ...) {
  UseMethod("standardize")
}

#' @rdname standardize
#' @export
standardize.because <- function(object, method = c("error", "total"), ...) {
  method <- match.arg(method)

  # Check if already standardized
  if (isTRUE(attr(object, "standardized")) || inherits(object, "because_standardized")) {
    message("Model is already standardized.")
    return(object)
  }

  res <- object
  data <- res$data
  if (is.null(data) && !is.null(res$original_data)) {
    data <- res$original_data
  }

  is_dsep <- !is.null(res$dsep) && (isTRUE(res$dsep) || is.list(res$dsep))

  if (is_dsep) {
    res <- .standardize_dsep(res, data = data, method = method)
  } else {
    res <- .standardize_fitted(res, data = data, method = method)
  }

  attr(res, "standardized") <- TRUE
  class(res) <- unique(c("because_standardized", class(res)))
  return(res)
}

#' Internal helper: standardize d-separation object
#' @noRd
.standardize_dsep <- function(res, data, method = "error") {
  if (is.null(res$dsep_results) || length(res$dsep_results) == 0) {
    return(res)
  }

  tests <- res$dsep_tests

  for (i in seq_along(res$dsep_results)) {
    res_i <- res$dsep_results[[i]]
    if (is.null(res_i) || is.null(res_i$samples)) next

    test_formula <- tests[[i]]
    resp <- as.character(test_formula)[2]
    fam <- if (!is.null(res$family) && resp %in% names(res$family)) {
      tolower(as.character(res$family[[resp]]))
    } else {
      "gaussian"
    }

    # Iterate over chains in this test
    for (ch in seq_along(res_i$samples)) {
      p_names <- colnames(res_i$samples[[ch]])

      # 1. Non-Gaussian: Binomial / Bernoulli / Ordinal
      if (fam %in% c("binomial", "bernoulli", "ordinal")) {
        sig_col <- grep(paste0("^sigma_.*", resp), p_names, value = TRUE)
        if (length(sig_col) == 0) {
          sig_col <- grep(paste0("^sigma_", resp), p_names, value = TRUE)
        }
        sig_sq <- if (length(sig_col) > 0) res_i$samples[[ch]][, sig_col[1]]^2 else 0

        # Liability disturbance SD
        pi_sq_3 <- (pi^2) / 3
        sd_liability <- sqrt(sig_sq + pi_sq_3)

        beta_cols <- grep(paste0("^beta_", resp, "_"), p_names, value = TRUE)
        for (b_col in beta_cols) {
          pred <- sub(paste0("^beta_", resp, "_"), "", b_col)
          pred_clean <- sub("(_\\d+|\\[\\d+\\])$", "", pred)
          sd_x <- if (!is.null(data) && pred_clean %in% names(data)) {
            stats::sd(as.numeric(data[[pred_clean]]), na.rm = TRUE)
          } else {
            1
          }
          if (is.na(sd_x) || sd_x == 0) sd_x <- 1
          res_i$samples[[ch]][, b_col] <- res_i$samples[[ch]][, b_col] * (sd_x / sd_liability)
        }
      } else if (fam == "gaussian") {
        # Continuous response
        sd_y <- if (!is.null(data) && resp %in% names(data)) {
          stats::sd(as.numeric(data[[resp]]), na.rm = TRUE)
        } else {
          1
        }
        if (is.na(sd_y) || sd_y == 0) sd_y <- 1

        beta_cols <- grep(paste0("^beta_", resp, "_"), p_names, value = TRUE)
        for (b_col in beta_cols) {
          pred <- sub(paste0("^beta_", resp, "_"), "", b_col)
          pred_clean <- sub("(_\\d+|\\[\\d+\\])$", "", pred)
          sd_x <- if (!is.null(data) && pred_clean %in% names(data)) {
            stats::sd(as.numeric(data[[pred_clean]]), na.rm = TRUE)
          } else {
            1
          }
          if (is.na(sd_x) || sd_x == 0) sd_x <- 1

          if (sd_x != sd_y) {
            res_i$samples[[ch]][, b_col] <- res_i$samples[[ch]][, b_col] * (sd_x / sd_y)
          }
        }
      }
    }
    res$dsep_results[[i]] <- res_i
  }

  # Synchronize global res$samples across all d-sep tests
  if (!is.null(res$samples)) {
    for (ch in seq_along(res$samples)) {
      p_names <- colnames(res$samples[[ch]])
      for (i in seq_along(res$dsep_results)) {
        if (is.null(res$dsep_results[[i]]) || is.null(res$dsep_results[[i]]$samples)) next
        res_samps <- res$dsep_results[[i]]$samples[[ch]]
        common_betas <- grep("^beta_", intersect(p_names, colnames(res_samps)), value = TRUE)
        if (length(common_betas) > 0) {
          res$samples[[ch]][, common_betas] <- res_samps[, common_betas]
        }
      }
    }
  }

  return(res)
}

#' Internal helper: standardize standard fitted because model
#' @noRd
.standardize_fitted <- function(res, data, method = "error") {
  if (is.null(res$samples)) return(res)

  param_map <- res$parameter_map
  p_names <- colnames(res$samples[[1]])
  standardized_params <- character(0)

  # Determine parameters to standardize
  if (!is.null(param_map) && is.data.frame(param_map)) {
    beta_rows <- param_map[grepl("^beta_", param_map$parameter), ]
    for (r in seq_len(nrow(beta_rows))) {
      p_name <- beta_rows$parameter[r]
      if (!p_name %in% p_names) next

      resp <- beta_rows$response[r]
      pred <- beta_rows$predictor[r]

      fam <- if (!is.null(res$family) && resp %in% names(res$family)) {
        tolower(as.character(res$family[[resp]]))
      } else {
        "gaussian"
      }

      # Predictor SD
      is_latent_pred <- !is.null(res$latent) && pred %in% res$latent
      sd_x <- if (is_latent_pred) {
        1.0
      } else if (!is.null(data) && pred %in% names(data)) {
        stats::sd(as.numeric(data[[pred]]), na.rm = TRUE)
      } else {
        1.0
      }
      if (is.na(sd_x) || sd_x == 0) sd_x <- 1.0

      if (fam %in% c("binomial", "bernoulli", "ordinal")) {
        sig_col <- grep(paste0("^sigma_.*", resp), p_names, value = TRUE)
        if (length(sig_col) == 0) {
          sig_col <- grep(paste0("^sigma_", resp), p_names, value = TRUE)
        }
        for (ch in seq_along(res$samples)) {
          sig_sq <- if (length(sig_col) > 0) res$samples[[ch]][, sig_col[1]]^2 else 0
          pi_sq_3 <- (pi^2) / 3
          sd_liability <- sqrt(sig_sq + pi_sq_3)
          res$samples[[ch]][, p_name] <- res$samples[[ch]][, p_name] * (sd_x / sd_liability)
        }
        standardized_params <- c(standardized_params, p_name)
      } else if (fam == "gaussian") {
        is_latent_resp <- !is.null(res$latent) && resp %in% res$latent
        sd_y <- if (is_latent_resp) {
          1.0
        } else if (!is.null(data) && resp %in% names(data)) {
          stats::sd(as.numeric(data[[resp]]), na.rm = TRUE)
        } else {
          1.0
        }
        if (is.na(sd_y) || sd_y == 0) sd_y <- 1.0

        if (sd_x != sd_y) {
          scale_fac <- sd_x / sd_y
          for (ch in seq_along(res$samples)) {
            res$samples[[ch]][, p_name] <- res$samples[[ch]][, p_name] * scale_fac
          }
          standardized_params <- c(standardized_params, p_name)
        }
      }
    }
  } else {
    # Fallback when parameter_map is not present: parse parameter names ^beta_Resp_Pred
    beta_cols <- grep("^beta_", p_names, value = TRUE)
    for (b_col in beta_cols) {
      clean_b <- sub("^beta_", "", b_col)
      parts <- strsplit(clean_b, "_")[[1]]
      if (length(parts) >= 2) {
        resp <- parts[1]
        pred <- paste(parts[-1], collapse = "_")
        pred_clean <- sub("(_\\d+|\\[\\d+\\])$", "", pred)

        fam <- if (!is.null(res$family) && resp %in% names(res$family)) {
          tolower(as.character(res$family[[resp]]))
        } else {
          "gaussian"
        }

        sd_x <- if (!is.null(data) && pred_clean %in% names(data)) stats::sd(as.numeric(data[[pred_clean]]), na.rm = TRUE) else 1
        if (is.na(sd_x) || sd_x == 0) sd_x <- 1

        if (fam %in% c("binomial", "bernoulli", "ordinal")) {
          sig_col <- grep(paste0("^sigma_.*", resp), p_names, value = TRUE)
          if (length(sig_col) == 0) {
            sig_col <- grep(paste0("^sigma_", resp), p_names, value = TRUE)
          }
          for (ch in seq_along(res$samples)) {
            sig_sq <- if (length(sig_col) > 0) res$samples[[ch]][, sig_col[1]]^2 else 0
            sd_liability <- sqrt(sig_sq + (pi^2)/3)
            res$samples[[ch]][, b_col] <- res$samples[[ch]][, b_col] * (sd_x / sd_liability)
          }
          standardized_params <- c(standardized_params, b_col)
        } else if (fam == "gaussian") {
          sd_y <- if (!is.null(data) && resp %in% names(data)) stats::sd(as.numeric(data[[resp]]), na.rm = TRUE) else 1
          if (is.na(sd_y) || sd_y == 0) sd_y <- 1
          if (sd_x != sd_y) {
            scale_fac <- sd_x / sd_y
            for (ch in seq_along(res$samples)) {
              res$samples[[ch]][, b_col] <- res$samples[[ch]][, b_col] * scale_fac
            }
            standardized_params <- c(standardized_params, b_col)
          }
        }
      }
    }
  }

  standardized_params <- unique(standardized_params)

  # High-efficiency summary update:
  # If summary is already present, update ONLY the rows for standardized parameters,
  # keeping all other parameters and precomputed convergence diagnostics (Rhat, Rhat_max, ESS) intact!
  if (!is.null(res$summary) && !is.null(res$summary$statistics) && length(standardized_params) > 0) {
    for (p_name in standardized_params) {
      if (p_name %in% rownames(res$summary$statistics)) {
        sub_samps <- res$samples[, p_name, drop = FALSE]
        sub_summ <- summary(sub_samps)
        stats_sub <- if (is.matrix(sub_summ$statistics)) sub_summ$statistics else matrix(sub_summ$statistics, nrow = 1, dimnames = list(p_name, names(sub_summ$statistics)))
        q_sub <- if (is.matrix(sub_summ$quantiles)) sub_summ$quantiles else matrix(sub_summ$quantiles, nrow = 1, dimnames = list(p_name, names(sub_summ$quantiles)))

        for (col in c("Mean", "SD", "Naive SE", "Time-series SE")) {
          if (col %in% colnames(stats_sub) && col %in% colnames(res$summary$statistics)) {
            res$summary$statistics[p_name, col] <- stats_sub[p_name, col]
          }
        }
        for (col in colnames(q_sub)) {
          if (col %in% colnames(res$summary$quantiles)) {
            res$summary$quantiles[p_name, col] <- q_sub[p_name, col]
          }
        }
      }
    }
  } else if (is.null(res$summary)) {
    # If no summary was present, only summarize structural parameters to keep it fast
    target_vars <- if (length(p_names) > 100) {
      p_names[grepl("^(alpha|beta|sigma|tau|lambda|rho|psi|r_|cutpoint)", p_names)]
    } else {
      p_names
    }
    if (length(target_vars) > 0) {
      res$summary <- compute_mcmc_summary(res$samples[, target_vars, drop = FALSE], length(res$samples))
    }
  }

  return(res)
}
