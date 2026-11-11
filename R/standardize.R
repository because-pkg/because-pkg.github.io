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
#'     \deqn{\beta^* = \frac{\beta \times \text{SD}(X)}{\sqrt{\sum \sigma_{\text{random}}^2 + \pi^2/3}}}
#'     where \eqn{\pi^2/3 \approx 3.289868} is the theoretical disturbance variance of the standard logistic link.
#'   \item \strong{Count (Poisson / Negative Binomial)}: Raw coefficients sit on the log-rate scale
#'     (\eqn{\ln \lambda}). Standardized onto the latent log / rate scale via:
#'     \itemize{
#'       \item \code{"error"} (default): \eqn{\beta^* = \frac{\beta \times \text{SD}(X)}{\sqrt{\sum \sigma_{\text{random}}^2 + \ln(1 + 1/\bar{Y})}}}
#'         following the GLMM latent distribution variance formulation of Nakagawa et al. (2017).
#'       \item \code{"total"}: Incorporates the empirical variance of the fixed linear predictor
#'         \eqn{\text{Var}(\hat{\eta}_{\text{fixed}})} alongside random effects and disturbance variance.
#'       \item \code{"mplus"}: Standardizes by the standard deviation of the linear predictions
#'         \eqn{s_{xb} = \sqrt{\text{Var}(\hat{\eta}_{\text{fixed}}) + \sum \sigma_{\text{random}}^2}},
#'         matching Mplus \code{StdYX} for Poisson regression (UCLA OARC FAQ).
#'     }
#'   \item \strong{Binary Predictors}: Slopes pointing from a binary indicator (e.g. \code{Migrant} \eqn{\in \{0, 1\}})
#'     to a response are kept on their natural 0-to-1 step when \code{binary_predictors = "raw"} (default),
#'     placing the effect on the standardized outcome scale (directly analogous to Cohen's \eqn{d}; Gelman 2008,
#'     Grace & Bollen 2005). When \code{binary_predictors = "sd"}, slopes are multiplied by \eqn{s_X = \sqrt{p(1 - p)}} (Menard 2004).
#'   \item \strong{d-separation Objects}: Evaluates and rescales each auxiliary independence claim
#'     across all d-sep sub-models, enabling a uniform \code{rope = c(-0.10, 0.10)} benchmark across all tests.
#' }
#'
#' @param object A fitted `because` model object or a `because` d-separation test object.
#' @param method Character; liability/outcome variance calculation method:
#'   \itemize{
#'     \item \code{"error"} (default): Disturbance error + structural/random effects variance
#'       (\eqn{\sum \sigma_{\text{random}}^2 + \pi^2/3} for binomial;
#'        \eqn{\sum \sigma_{\text{random}}^2 + \ln(1 + 1/\bar{Y})} for Poisson; Nakagawa et al. 2017).
#'     \item \code{"total"}: Full McKelvey-Zavoina / Nakagawa total variance incorporating the empirical variance
#'       of the fixed linear predictor \eqn{\text{Var}(\hat{\eta}_{\text{fixed}})}.
#'     \item \code{"mplus"}: Standardizes by the standard deviation of the linear predictions
#'       \eqn{s_{xb} = \sqrt{\text{Var}(\hat{\eta}_{\text{fixed}}) + \sum \sigma_{\text{random}}^2}},
#'       matching Mplus \code{StdYX} for Poisson regression (UCLA OARC FAQ).
#'   }
#' @param binary_predictors Character; scaling for binary/dummy predictors: `"raw"` (default;
#'   semi-standardization, keeping binary predictors on the 0-to-1 step so that the standardized slope
#'   represents Cohen's \eqn{d} on the standardized outcome scale; Gelman 2008, Bollen 1989, Grace & Bollen 2005)
#'   or `"sd"` (full standardization, scaling by \eqn{s_X = \sqrt{p(1-p)}}; Menard 2004).
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
#' Gelman, A. (2008). Scaling regression inputs by dividing by two standard deviations.
#' \emph{Statistics in Medicine}, 27(15), 2865–2873.
#'
#' Grace, J. B., & Bollen, K. A. (2005). Interpreting the results from multiple regression and structural equation models.
#' \emph{Ecology}, 86(11), 283–295.
#'
#' Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The coefficient of determination \eqn{R^2}
#' and intra-class correlation coefficient from generalized linear mixed-effects models revisited and expanded.
#' \emph{Journal of the Royal Society Interface}, 14(134), 20170213.
#'
#' UCLA: Statistical Consulting Group. How does Mplus calculate the standardized coefficients based on a Poisson model?
#' \url{https://stats.oarc.ucla.edu/mplus/faq/how-does-mplus-calculate-the-standardized-coefficients-based-on-a-poisson-model/}
#'
#' @export
standardize <- function(object, ...) {
  UseMethod("standardize")
}

#' @rdname standardize
#' @export
standardize.because <- function(object, method = c("error", "total", "mplus"), binary_predictors = c("raw", "sd"), ...) {
  method <- match.arg(method)
  binary_predictors <- match.arg(binary_predictors)

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
    res <- .standardize_dsep(res, data = data, method = method, binary_predictors = binary_predictors)
  } else {
    res <- .standardize_fitted(res, data = data, method = method, binary_predictors = binary_predictors)
  }

  attr(res, "standardized") <- TRUE
  attr(res, "binary_predictors") <- binary_predictors
  class(res) <- unique(c("because_standardized", class(res)))
  return(res)
}

#' Helper to check if a variable is binary (0/1 or 2-level factor)
#' @noRd
.is_binary_var <- function(var_name, data) {
  if (is.null(data) || !var_name %in% names(data)) return(FALSE)
  x <- stats::na.omit(data[[var_name]])
  if (is.factor(x) || is.logical(x)) {
    return(length(unique(x)) <= 2)
  }
  if (is.numeric(x)) {
    u_vals <- unique(x)
    return(length(u_vals) <= 2 && all(u_vals %in% c(0, 1)))
  }
  FALSE
}

#' Helper to determine variable family with fallbacks
#' @noRd
.get_var_family <- function(res, resp, data = NULL) {
  fam <- if (!is.null(res$family) && resp %in% names(res$family)) {
    tolower(as.character(res$family[[resp]]))
  } else if (!is.null(res$input$family) && resp %in% names(res$input$family)) {
    tolower(as.character(res$input$family[[resp]]))
  } else if (!is.null(res$py_result$family) && resp %in% names(res$py_result$family)) {
    tolower(as.character(res$py_result$family[[resp]]))
  } else {
    NULL
  }

  if (!is.null(fam) && length(fam) > 0 && fam != "") {
    return(fam)
  }

  # Inspect data if available (e.g. for pre-existing objects loaded from disk)
  if (!is.null(data) && resp %in% names(data)) {
    y <- stats::na.omit(data[[resp]])
    if (is.factor(y) || is.logical(y)) {
      if (is.ordered(y)) return("ordinal")
      return("binomial")
    }
    if (is.numeric(y)) {
      u_vals <- unique(y)
      if (length(u_vals) > 0 && all(u_vals %in% c(0, 1))) {
        return("binomial")
      }
      if (length(y) > 0 && all(y >= 0, na.rm = TRUE) && all(y == round(y), na.rm = TRUE) && length(u_vals) > 2) {
        return("poisson")
      }
    }
  }

  return("gaussian")
}

#' Helper to calculate variance of fixed linear predictor (Var(X * beta))
#' @noRd
.get_fixed_variance <- function(samples_ch, resp, p_names, param_map = NULL, data = NULL) {
  if (is.null(data)) return(0)

  pred_names <- character(0)
  beta_cols <- character(0)
  if (!is.null(param_map) && is.data.frame(param_map)) {
    matching_pm <- param_map[param_map$response == resp & grepl("^beta_", param_map$parameter), ]
    if (nrow(matching_pm) > 0) {
      for (r in seq_len(nrow(matching_pm))) {
        p_col <- matching_pm$parameter[r]
        pred_var <- matching_pm$predictor[r]
        if (p_col %in% p_names && pred_var %in% names(data)) {
          pred_names <- c(pred_names, pred_var)
          beta_cols <- c(beta_cols, p_col)
        }
      }
    }
  }

  if (length(beta_cols) == 0) {
    b_cand <- grep(paste0("^beta_", resp, "_"), p_names, value = TRUE)
    for (bc in b_cand) {
      pred_cand <- sub(paste0("^beta_", resp, "_"), "", bc)
      pred_clean <- sub("(_\\d+|\\[\\d+\\])$", "", pred_cand)
      if (pred_clean %in% names(data)) {
        pred_names <- c(pred_names, pred_clean)
        beta_cols <- c(beta_cols, bc)
      }
    }
  }

  if (length(beta_cols) == 0) return(0)

  u_idx <- !duplicated(beta_cols)
  beta_cols <- beta_cols[u_idx]
  pred_names <- pred_names[u_idx]

  x_mat <- as.matrix(data[, pred_names, drop = FALSE])
  mode(x_mat) <- "numeric"
  cov_x <- stats::cov(x_mat, use = "pairwise.complete.obs")
  cov_x[is.na(cov_x)] <- 0

  beta_mat <- as.matrix(samples_ch[, beta_cols, drop = FALSE])
  var_fixed <- rowSums((beta_mat %*% cov_x) * beta_mat)
  var_fixed[is.na(var_fixed) | var_fixed < 0] <- 0
  return(var_fixed)
}

#' Helper to compute liability scale standard deviation for non-Gaussian response
#' @noRd
.get_liability_sd <- function(samples_ch, resp, p_names, param_map = NULL,
                              family = "binomial", data = NULL, method = "error") {
  fam <- tolower(as.character(family))

  # Distribution disturbance variance on the latent/link scale
  if (fam %in% c("poisson", "negbinomial")) {
    if (method == "mplus") {
      dist_var <- 0 # Pure Mplus StdYX on linear predictor without distribution noise
    } else {
      mean_y <- if (!is.null(data) && resp %in% names(data)) {
        mean(as.numeric(data[[resp]]), na.rm = TRUE)
      } else {
        1.0
      }
      if (is.na(mean_y) || mean_y <= 0) mean_y <- 1.0
      dist_var <- log(1 + 1 / mean_y)
    }
  } else {
    # Binomial / Bernoulli / Ordinal (logistic theoretical disturbance pi^2/3)
    dist_var <- (pi^2) / 3
  }

  # Structural / Random effects variance
  sig_cols <- character(0)
  if (!is.null(param_map) && is.data.frame(param_map)) {
    matching_pm <- param_map[param_map$response == resp & grepl("^sigma_", param_map$parameter), "parameter"]
    sig_cols <- intersect(matching_pm, p_names)
    sig_cols <- sig_cols[!grepl("^sigma_tot_", sig_cols)]
  }

  if (length(sig_cols) == 0) {
    # Match structural sigmas for this response (e.g. sigma_Migrant_multiPhylo, sigma_Migrant)
    sig_cols <- grep(paste0("^sigma_(.*_)?", resp, "(_.*)?$"), p_names, value = TRUE)
    sig_cols <- sig_cols[!grepl("^sigma_tot_", sig_cols)]
  }

  struct_var <- 0
  if (length(sig_cols) > 0) {
    for (sc in sig_cols) {
      struct_var <- struct_var + samples_ch[, sc]^2
    }
  }

  # Fixed predictor variance (for method = "total" or "mplus")
  fixed_var <- 0
  if (method %in% c("total", "mplus")) {
    fixed_var <- .get_fixed_variance(samples_ch, resp, p_names, param_map, data)
  }

  tot_var <- fixed_var + struct_var + dist_var
  tot_var <- pmax(tot_var, 1e-8)
  return(sqrt(tot_var))
}

#' Internal helper: standardize d-separation object
#' @noRd
.standardize_dsep <- function(res, data, method = "error", binary_predictors = "raw") {
  if (is.null(res$dsep_results) || length(res$dsep_results) == 0) {
    return(res)
  }

  tests <- res$dsep_tests

  for (i in seq_along(res$dsep_results)) {
    res_i <- res$dsep_results[[i]]
    if (is.null(res_i) || is.null(res_i$samples)) next

    test_formula <- tests[[i]]
    resp <- as.character(test_formula)[2]
    fam <- .get_var_family(res, resp, data)

    # Iterate over chains in this test
    for (ch in seq_along(res_i$samples)) {
      p_names <- colnames(res_i$samples[[ch]])

      # 1. Non-Gaussian: Binomial / Bernoulli / Ordinal / Poisson / NegBinomial
      if (fam %in% c("binomial", "bernoulli", "ordinal", "poisson", "negbinomial")) {
        sd_liability <- .get_liability_sd(res_i$samples[[ch]], resp, p_names, res_i$parameter_map,
                                          family = fam, data = data, method = method)

        beta_cols <- grep(paste0("^beta_", resp, "_"), p_names, value = TRUE)
        for (b_col in beta_cols) {
          pred <- sub(paste0("^beta_", resp, "_"), "", b_col)
          pred_clean <- sub("(_\\d+|\\[\\d+\\])$", "", pred)
          is_binary_pred <- .is_binary_var(pred_clean, data)
          sd_x <- if (is_binary_pred && binary_predictors == "raw") {
            1.0
          } else if (!is.null(data) && pred_clean %in% names(data)) {
            stats::sd(as.numeric(data[[pred_clean]]), na.rm = TRUE)
          } else {
            1.0
          }
          if (is.na(sd_x) || sd_x == 0) sd_x <- 1.0
          res_i$samples[[ch]][, b_col] <- res_i$samples[[ch]][, b_col] * (sd_x / sd_liability)
        }
      } else if (fam == "gaussian") {
        # Continuous response
        sd_y <- if (!is.null(data) && resp %in% names(data)) {
          stats::sd(as.numeric(data[[resp]]), na.rm = TRUE)
        } else {
          1.0
        }
        if (is.na(sd_y) || sd_y == 0) sd_y <- 1.0

        beta_cols <- grep(paste0("^beta_", resp, "_"), p_names, value = TRUE)
        for (b_col in beta_cols) {
          pred <- sub(paste0("^beta_", resp, "_"), "", b_col)
          pred_clean <- sub("(_\\d+|\\[\\d+\\])$", "", pred)
          is_binary_pred <- .is_binary_var(pred_clean, data)
          sd_x <- if (is_binary_pred && binary_predictors == "raw") {
            1.0
          } else if (!is.null(data) && pred_clean %in% names(data)) {
            stats::sd(as.numeric(data[[pred_clean]]), na.rm = TRUE)
          } else {
            1.0
          }
          if (is.na(sd_x) || sd_x == 0) sd_x <- 1.0

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
.standardize_fitted <- function(res, data, method = "error", binary_predictors = "raw") {
  if (is.null(res$samples)) return(res)

  # Keep an unstandardized copy of samples so liability SD calculations are not corrupted by in-place updates
  orig_samples <- res$samples
  param_map <- res$parameter_map
  p_names <- colnames(res$samples[[1]])
  standardized_params <- character(0)

  # Precompute liability SD per response and chain to ensure efficiency and consistency across predictors
  liability_sd_cache <- list()

  # Determine parameters to standardize
  if (!is.null(param_map) && is.data.frame(param_map)) {
    beta_rows <- param_map[grepl("^beta_", param_map$parameter), ]
    for (r in seq_len(nrow(beta_rows))) {
      p_name <- beta_rows$parameter[r]
      if (!p_name %in% p_names) next

      resp <- beta_rows$response[r]
      pred <- beta_rows$predictor[r]

      fam <- .get_var_family(res, resp, data)

      # Predictor SD
      is_latent_pred <- !is.null(res$latent) && pred %in% res$latent
      is_binary_pred <- .is_binary_var(pred, data)
      sd_x <- if (is_latent_pred) {
        1.0
      } else if (is_binary_pred && binary_predictors == "raw") {
        1.0 # Keep on 0-to-1 step for Cohen's d (Gelman 2008, Grace & Bollen 2005)
      } else if (!is.null(data) && pred %in% names(data)) {
        stats::sd(as.numeric(data[[pred]]), na.rm = TRUE)
      } else {
        1.0
      }
      if (is.na(sd_x) || sd_x == 0) sd_x <- 1.0

      if (fam %in% c("binomial", "bernoulli", "ordinal", "poisson", "negbinomial")) {
        for (ch in seq_along(res$samples)) {
          cache_key <- paste0(resp, "_ch", ch)
          if (is.null(liability_sd_cache[[cache_key]])) {
            liability_sd_cache[[cache_key]] <- .get_liability_sd(
              orig_samples[[ch]], resp, p_names, res$parameter_map,
              family = fam, data = data, method = method
            )
          }
          sd_liability <- liability_sd_cache[[cache_key]]
          res$samples[[ch]][, p_name] <- orig_samples[[ch]][, p_name] * (sd_x / sd_liability)
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

        fam <- .get_var_family(res, resp, data)

        is_binary_pred <- .is_binary_var(pred_clean, data)
        sd_x <- if (is_binary_pred && binary_predictors == "raw") {
          1.0 # Keep on 0-to-1 step for Cohen's d
        } else if (!is.null(data) && pred_clean %in% names(data)) {
          stats::sd(as.numeric(data[[pred_clean]]), na.rm = TRUE)
        } else {
          1.0
        }
        if (is.na(sd_x) || sd_x == 0) sd_x <- 1.0

        if (fam %in% c("binomial", "bernoulli", "ordinal", "poisson", "negbinomial")) {
          for (ch in seq_along(res$samples)) {
            cache_key <- paste0(resp, "_ch", ch)
            if (is.null(liability_sd_cache[[cache_key]])) {
              liability_sd_cache[[cache_key]] <- .get_liability_sd(
                orig_samples[[ch]], resp, p_names, NULL,
                family = fam, data = data, method = method
              )
            }
            sd_liability <- liability_sd_cache[[cache_key]]
            res$samples[[ch]][, b_col] <- orig_samples[[ch]][, b_col] * (sd_x / sd_liability)
          }
          standardized_params <- c(standardized_params, b_col)
        } else if (fam == "gaussian") {
          sd_y <- if (!is.null(data) && resp %in% names(data)) stats::sd(as.numeric(data[[resp]]), na.rm = TRUE) else 1.0
          if (is.na(sd_y) || sd_y == 0) sd_y <- 1.0
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
