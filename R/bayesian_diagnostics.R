#' @title Purely Bayesian Assumption Diagnostics for because Models
#' @name bayesian_diagnostics
#' @description
#' A comprehensive, purely Bayesian suite of diagnostic functions for testing
#' regression and DAG/SEM assumptions in \code{because} models.
#'
#' In strict accordance with Bayesian data analysis philosophy (Gelman et al. 2013;
#' McElreath 2020), this suite excludes all frequentist null hypothesis testing
#' and \eqn{p}-values. Instead, diagnostics evaluate posterior parameter covariance
#' geometry and posterior predictive discrepancies.
#'
#' Included diagnostics:
#' \itemize{
#'   \item \code{\link{check_collinearity}} / \code{\link{bvif}}: Leamer's Bayesian Variance
#'         Inflation Factor and posterior parameter correlation matrix.
#'   \item \code{\link{check_normality}} / \code{\link{plot_qq}}: Bayesian Q-Q check with
#'         empirical posterior predictive simulation envelopes and skewness/kurtosis discrepancies.
#'   \item \code{\link{check_homoscedasticity}}: Binned posterior residual standard deviations
#'         across fitted values and variance ratio discrepancies.
#'   \item \code{\link{check_residual_structure}}: Posterior distribution of the Variance Partition
#'         Coefficient (\eqn{\lambda_{\text{Bayes}}}) for phylogenetic or spatial random effects.
#'   \item \code{\link{check_model}}: Unified multi-panel Bayesian diagnostic dashboard.
#' }
#'
#' @references
#' Leamer, E. E. (1973). Multicollinearity in regression analysis: An alternative view.
#'   \emph{Review of Economics and Statistics}, 55(3), 371-380.
#'
#' Leamer, E. E. (1978). \emph{Specification Searches: Ad Hoc Inference with Nonexperimental Data}.
#'   John Wiley & Sons.
#'
#' McElreath, R. (2020). \emph{Statistical Rethinking: A Bayesian Course with Examples in R and Stan} (2nd ed.).
#'   CRC Press.
#'
#' Gelman, A., Meng, X.-L., & Stern, H. (1996). Posterior predictive assessment of model fitness
#'   via realized discrepancies. \emph{Statistica Sinica}, 6(4), 733-760.
#'
#' Landwehr, J. M., Pregibon, D., & Shoemaker, A. C. (1984). Graphical methods for assessing logistic regression models.
#'   \emph{Journal of the American Statistical Association}, 79(385), 61-71.
#'
#' Hadfield, J. D. (2010). MCMC methods for multi-response generalized linear mixed models: The MCMCglmm R package.
#'   \emph{Journal of Statistical Software}, 33(2), 1-22.
#'
#' @importFrom stats cov cor solve qnorm quantile median sd
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon geom_hline geom_col geom_errorbar labs theme_minimal coord_flip facet_wrap
NULL

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==============================================================================
# 1. Collinearity & Bayesian Variance Inflation (BVIF)
# ==============================================================================

#' @title Bayesian Collinearity & Variance Inflation Diagnostics
#' @rdname check_collinearity
#' @description
#' Evaluates collinearity among predictors in each equation of a fitted \code{because} model.
#'
#' Computes the **Posterior Parameter Correlation Matrix** \eqn{\text{Cor}(\boldsymbol{\beta} \mid y)}
#' and **Leamer's Bayesian Variance Inflation Factor (BVIF)**:
#' \deqn{\text{BVIF}_j = [\boldsymbol{\Sigma}_{\boldsymbol{\beta}}]_{jj} \cdot [\boldsymbol{\Sigma}_{\boldsymbol{\beta}}^{-1}]_{jj} = \frac{1}{1 - R^2_{\beta_j \mid \boldsymbol{\beta}_{-j}}}}
#' where \eqn{\boldsymbol{\Sigma}_{\boldsymbol{\beta}}} is the posterior covariance matrix of the regression parameters.
#'
#' The square root of BVIF gives the **SE Multiplier**, quantifying the expansion in posterior
#' uncertainty attributable to shared informational content between predictors.
#'
#' @param object A \code{because} model fit object.
#' @param equation Optional character string or integer specifying a single equation to evaluate.
#'   If \code{NULL} (default), evaluates all equations with 2 or more predictors.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"because_collinearity"}, containing summary tables and correlation matrices.
#'
#' @details
#' \strong{Causal Notice:} In a causal DAG, elevated variance inflation between confounders or mediators
#' is an expected consequence of network topology. Variables required by the DAG's adjustment set
#' should \strong{not} be dropped based on BVIF, as doing so introduces confounding bias.
#'
#' @export
check_collinearity.because <- function(object, equation = NULL, ...) {
  if (!inherits(object, "because")) {
    stop("object must be of class 'because'")
  }
  
  samples_mat <- as.matrix(object$samples)
  pm <- object$parameter_map
  eqs <- object$equations
  
  results_list <- list()
  cor_matrices <- list()
  
  for (j in seq_along(eqs)) {
    eq <- eqs[[j]]
    eq_form <- if (is.list(eq) && "formula" %in% names(eq)) eq$formula else stats::formula(eq)
    resp <- as.character(eq_form)[2]
    
    if (!is.null(equation)) {
      if (is.character(equation) && equation != resp && !grepl(equation, deparse(eq_form))) next
      if (is.numeric(equation) && equation != j) next
    }
    
    # Identify fixed structural predictors (excluding intercept and random effects)
    pm_sub <- pm[pm$response == resp & pm$type == "coefficient" & pm$predictor != "(Intercept)", ]
    
    # Filter to parameters present in samples
    avail_params <- pm_sub$parameter[pm_sub$parameter %in% colnames(samples_mat)]
    avail_preds  <- pm_sub$predictor[pm_sub$parameter %in% colnames(samples_mat)]
    
    if (length(avail_params) < 2) {
      next
    }
    
    # Extract posterior draws
    beta_mat <- samples_mat[, avail_params, drop = FALSE]
    colnames(beta_mat) <- avail_preds
    
    # Posterior covariance and correlation
    sigma_beta <- stats::cov(beta_mat)
    cor_beta   <- stats::cor(beta_mat)
    cor_matrices[[resp]] <- cor_beta
    
    # Inverse covariance (precision matrix)
    omega_beta <- tryCatch(
      stats::solve(sigma_beta),
      error = function(e) {
        # Fallback with slight regularization if singular
        stats::solve(sigma_beta + diag(1e-6, ncol(sigma_beta)))
      }
    )
    
    # Leamer BVIF: diag(Sigma) * diag(Omega)
    bvif_vals <- diag(sigma_beta) * diag(omega_beta)
    se_mult   <- sqrt(bvif_vals)
    
    # Identify strongest correlation partner
    max_cor_with <- character(length(avail_preds))
    max_cor_r    <- numeric(length(avail_preds))
    
    for (k in seq_along(avail_preds)) {
      cors_k <- cor_beta[k, -k]
      if (length(cors_k) > 0) {
        abs_cors <- abs(cors_k)
        max_idx  <- which.max(abs_cors)
        max_cor_with[k] <- names(cors_k)[max_idx]
        max_cor_r[k]    <- cors_k[max_idx]
      } else {
        max_cor_with[k] <- NA_character_
        max_cor_r[k]    <- 0
      }
    }
    
    df_resp <- data.frame(
      Response      = resp,
      Term          = avail_preds,
      BVIF          = round(bvif_vals, 2),
      SE_Multiplier = round(se_mult, 2),
      Max_Cor_With  = max_cor_with,
      Max_Cor_r     = round(max_cor_r, 2),
      stringsAsFactors = FALSE
    )
    results_list[[resp]] <- df_resp
  }
  
  if (length(results_list) == 0) {
    message("No equations with >= 2 structural predictors found.")
    return(invisible(NULL))
  }
  
  summary_table <- do.call(rbind, results_list)
  rownames(summary_table) <- NULL
  
  out <- list(
    summary      = summary_table,
    cor_matrices = cor_matrices,
    equations    = eqs
  )
  class(out) <- c("because_collinearity", "list")
  return(out)
}

#' @rdname check_collinearity
#' @export
bvif <- function(object, ...) {
  check_collinearity(object, ...)
}

#' @export
print.because_collinearity <- function(x, ...) {
  cat("\n=== Bayesian Collinearity & Variance Inflation Diagnostics ===\n")
  cat("Framework: Leamer (1973, 1978) Posterior Variance Inflation\n\n")
  
  print(format(x$summary, justify = "left"), row.names = FALSE)
  
  cat("\nLegend:\n")
  cat("  BVIF          : Bayesian Variance Inflation Factor (1 = orthogonal)\n")
  cat("  SE_Multiplier : Posterior uncertainty inflation factor (sqrt(BVIF))\n")
  cat("  Max_Cor_With  : Strongest posterior correlation partner\n")
  cat("  Max_Cor_r     : Strongest posterior correlation value\n\n")
  
  cat("Causal Guidance (McElreath 2020):\n")
  cat("  In a causal DAG, elevated variance inflation between confounders or mediators\n")
  cat("  is expected. Do not drop variables required by the causal DAG's adjustment set;\n")
  cat("  doing so introduces confounding bias. BVIF reflects precision loss (SE expansion),\n")
  cat("  not model invalidity.\n\n")
  invisible(x)
}

#' @export
plot.because_collinearity <- function(x, ...) {
  df <- x$summary
  if (nrow(df) == 0) return(invisible(NULL))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(paste0(Response, ": ", Term), BVIF), y = BVIF)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8, width = 0.6) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen") +
    ggplot2::geom_hline(yintercept = 5, linetype = "dotted", color = "darkorange") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = "Bayesian Variance Inflation Factors (BVIF)",
      subtitle = "Dashed = 1.0 (Orthogonal); Dotted = 5.0 (Moderate Inflation)",
      x        = "Predictor Path",
      y        = "BVIF (Posterior Variance Inflation)"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}


# ==============================================================================
# 2. Normality & Posterior Predictive Q-Q Envelopes
# ==============================================================================

#' @title Bayesian Normality Check & Q-Q Envelope
#' @rdname check_normality
#' @description
#' Evaluates distributional and normality assumptions using **Randomized Quantile Residuals**
#' (Dunn & Smyth 1996; Hartig 2022) and **Posterior Predictive Credible Envelopes** (Landwehr et al. 1984).
#'
#' Generates an empirical simulation envelope (50%, 80%, and 95% bands) directly from posterior
#' predictive replications \eqn{y^{\text{rep}}}, and evaluates sample skewness and kurtosis
#' test quantities against their posterior predictive distributions.
#'
#' @param object A \code{because} model fit object.
#' @param resp Character string; name of the response variable. If \code{NULL}, takes the first response.
#' @param ndraws Integer; number of posterior predictive draws to simulate for the envelope. Defaults to 250.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"because_normality"} containing the ggplot object, envelope data,
#'   and skewness/kurtosis posterior tail probabilities.
#'
#' @export
check_normality.because <- function(object, resp = NULL, ndraws = 250, ...) {
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  # 1. Observed Quantile Residuals
  res_obs <- residuals(object, resp = resp, type = "quantile", ndraws = ndraws)
  n <- length(res_obs)
  
  # 2. Generate Posterior Predictive Replicates for Empirical Envelope
  yrep <- posterior_predict(object, resp = resp, ndraws = min(ndraws, 100))
  n_rep <- nrow(yrep)
  
  # Theoretical quantiles
  theo_q <- stats::qnorm((seq_len(n) - 0.375) / (n + 0.25))
  obs_sorted <- sort(res_obs)
  
  # Compute sorted replicates under the model
  rep_sorted_mat <- matrix(0, nrow = n_rep, ncol = n)
  for (b in seq_len(n_rep)) {
    # Generate quantile residuals for replicated data
    u_b <- numeric(n)
    for (i in seq_len(n)) {
      u_b[i] <- mean(yrep[, i] < yrep[b, i]) + stats::runif(1) * mean(yrep[, i] == yrep[b, i])
      u_b[i] <- min(max(u_b[i], 1 / (2 * n_rep)), 1 - 1 / (2 * n_rep))
    }
    rep_sorted_mat[b, ] <- sort(stats::qnorm(u_b))
  }
  
  # Empirical Envelope Quantiles
  lower95 <- apply(rep_sorted_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  lower80 <- apply(rep_sorted_mat, 2, stats::quantile, probs = 0.10,  na.rm = TRUE)
  median_q <- apply(rep_sorted_mat, 2, stats::median, na.rm = TRUE)
  upper80 <- apply(rep_sorted_mat, 2, stats::quantile, probs = 0.90,  na.rm = TRUE)
  upper95 <- apply(rep_sorted_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  
  envelope_df <- data.frame(
    Theoretical = theo_q,
    Observed    = obs_sorted,
    Lower95     = lower95,
    Lower80     = lower80,
    Median      = median_q,
    Upper80     = upper80,
    Upper95     = upper95
  )
  
  # 3. Discrepancy Statistics: Skewness and Kurtosis
  skewness_fn <- function(x) {
    m2 <- mean((x - mean(x))^2)
    if (m2 == 0) return(0)
    mean((x - mean(x))^3) / (m2^(3/2))
  }
  kurtosis_fn <- function(x) {
    m2 <- mean((x - mean(x))^2)
    if (m2 == 0) return(0)
    mean((x - mean(x))^4) / (m2^2)
  }
  
  t_skew_obs <- skewness_fn(res_obs)
  t_kurt_obs <- kurtosis_fn(res_obs)
  
  t_skew_rep <- apply(rep_sorted_mat, 1, skewness_fn)
  t_kurt_rep <- apply(rep_sorted_mat, 1, kurtosis_fn)
  
  tail_p_skew <- mean(t_skew_rep >= t_skew_obs)
  tail_p_kurt <- mean(t_kurt_rep >= t_kurt_obs)
  
  discrepancies <- data.frame(
    Metric        = c("Skewness", "Kurtosis"),
    Observed      = c(round(t_skew_obs, 2), round(t_kurt_obs, 2)),
    Rep_Median    = c(round(stats::median(t_skew_rep), 2), round(stats::median(t_kurt_rep), 2)),
    Rep_95_CI     = c(
      paste0("[", round(stats::quantile(t_skew_rep, 0.025), 2), ", ", round(stats::quantile(t_skew_rep, 0.975), 2), "]"),
      paste0("[", round(stats::quantile(t_kurt_rep, 0.025), 2), ", ", round(stats::quantile(t_kurt_rep, 0.975), 2), "]")
    ),
    Tail_Prob     = c(round(tail_p_skew, 3), round(tail_p_kurt, 3)),
    stringsAsFactors = FALSE
  )
  
  # 4. Construct Plot
  p <- ggplot2::ggplot(envelope_df, ggplot2::aes(x = Theoretical)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower95, ymax = Upper95), fill = "skyblue", alpha = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower80, ymax = Upper80), fill = "skyblue", alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = Theoretical), linetype = "dashed", color = "darkgray") +
    ggplot2::geom_point(ggplot2::aes(y = Observed), color = "darkblue", size = 1.3, alpha = 0.8) +
    ggplot2::labs(
      title    = paste("Bayesian Q-Q Plot with Credible Envelope:", resp),
      subtitle = "Bands: 50% & 95% Posterior Predictive Envelopes (Landwehr et al. 1984)",
      x        = "Theoretical Normal Quantiles",
      y        = "Randomized Quantile Residuals"
    ) +
    ggplot2::theme_minimal()
  
  out <- list(
    plot          = p,
    discrepancies = discrepancies,
    envelope      = envelope_df,
    response      = resp
  )
  class(out) <- c("because_normality", "list")
  return(out)
}

#' @rdname check_normality
#' @export
plot_qq.because <- function(object, ...) {
  res <- check_normality(object, ...)
  return(res$plot)
}

#' @export
print.because_normality <- function(x, ...) {
  cat("\n=== Bayesian Normality Check (Randomized Quantile Residuals) ===\n")
  cat(paste("Response Variable:", x$response, "\n\n"))
  cat("Posterior Predictive Discrepancies (Gelman et al. 1996):\n")
  print(format(x$discrepancies, justify = "left"), row.names = FALSE)
  cat("\nNote: Tail_Prob = P(T(y_rep) >= T(y)). Extreme values (< 0.025 or > 0.975)\n")
  cat("indicate that observed skewness or kurtosis is atypical under model replications.\n\n")
  invisible(x)
}


# ==============================================================================
# 3. Homoscedasticity & Binned Posterior Variances
# ==============================================================================

#' @title Bayesian Homoscedasticity & Residual Dispersion Diagnostics
#' @rdname check_homoscedasticity
#' @description
#' Assesses variance constancy across fitted predictions \eqn{\hat{\mu}} in a purely Bayesian manner
#' (Gelman & Hill 2007).
#'
#' Partitions fitted values into quantiles and computes the posterior median and 95% Credible Interval
#' of residual standard deviations in each bin. Also computes the posterior distribution of the
#' variance ratio discrepancy:
#' \deqn{T_{\text{ratio}}^{(s)} = \frac{\text{Var}(e_{\text{top\_bin}}^{(s)})}{\text{Var}(e_{\text{bottom\_bin}}^{(s)})}}
#'
#' @param object A \code{because} model fit object.
#' @param resp Character string; name of the response variable. If \code{NULL}, takes the first response.
#' @param n_bins Integer; number of bins to partition fitted values into. Defaults to 5 (quintiles).
#' @param ndraws Integer; number of posterior draws to use. Defaults to 250.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"because_homoscedasticity"} containing diagnostic plots and summary statistics.
#'
#' @export
check_homoscedasticity.because <- function(object, resp = NULL, n_bins = 5, ndraws = 250, ...) {
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  # Fitted values and residuals
  mu_hat <- fitted(object, resp = resp, summary = TRUE)
  res_mat <- residuals(object, resp = resp, type = "response", summary = FALSE, ndraws = ndraws)
  
  n_s <- nrow(res_mat)
  n_obs <- ncol(res_mat)
  
  if (length(mu_hat) != n_obs) {
    mu_hat <- mu_hat[seq_len(n_obs)]
  }
  
  # Bin fitted values into quantiles
  breaks <- stats::quantile(mu_hat, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  # Ensure unique breaks
  breaks <- unique(breaks)
  if (length(breaks) < 3) {
    breaks <- seq(min(mu_hat, na.rm = TRUE), max(mu_hat, na.rm = TRUE), length.out = n_bins + 1)
  }
  bins <- cut(mu_hat, breaks = breaks, include.lowest = TRUE)
  bin_levels <- levels(bins)
  
  # Compute residual SD per bin across all posterior draws
  bin_sd_mat <- matrix(NA, nrow = n_s, ncol = length(bin_levels))
  colnames(bin_sd_mat) <- bin_levels
  
  for (k in seq_along(bin_levels)) {
    idx_k <- which(bins == bin_levels[k])
    if (length(idx_k) >= 2) {
      bin_sd_mat[, k] <- apply(res_mat[, idx_k, drop = FALSE], 1, stats::sd, na.rm = TRUE)
    }
  }
  
  # Binned summary
  bin_summary <- data.frame(
    Bin       = bin_levels,
    Median_SD = round(apply(bin_sd_mat, 2, stats::median, na.rm = TRUE), 3),
    Lower95   = round(apply(bin_sd_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE), 3),
    Upper95   = round(apply(bin_sd_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE), 3)
  )
  
  # Variance Ratio Discrepancy (Top bin vs Bottom bin)
  top_bin_idx <- length(bin_levels)
  var_ratio_samples <- (bin_sd_mat[, top_bin_idx]^2) / (bin_sd_mat[, 1]^2)
  var_ratio_median  <- stats::median(var_ratio_samples, na.rm = TRUE)
  var_ratio_ci      <- stats::quantile(var_ratio_samples, probs = c(0.025, 0.975), na.rm = TRUE)
  
  overall_sd <- stats::median(apply(res_mat, 1, stats::sd, na.rm = TRUE))
  
  # Plot 1: Residuals vs Fitted
  plot_df_res <- data.frame(Fitted = mu_hat, Residual = colMeans(res_mat, na.rm = TRUE))
  p1 <- ggplot2::ggplot(plot_df_res, ggplot2::aes(x = Fitted, y = Residual)) +
    ggplot2::geom_point(color = "darkslategray", alpha = 0.7, size = 1.3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick") +
    ggplot2::labs(
      title = paste("Residuals vs. Fitted:", resp),
      x = "Posterior Expected Fitted Value",
      y = "Posterior Mean Residual"
    ) +
    ggplot2::theme_minimal()
  
  # Plot 2: Binned Residual Variance
  bin_summary$Bin_Num <- seq_len(nrow(bin_summary))
  p2 <- ggplot2::ggplot(bin_summary, ggplot2::aes(x = factor(Bin_Num), y = Median_SD)) +
    ggplot2::geom_point(color = "navy", size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower95, ymax = Upper95), width = 0.2, color = "navy") +
    ggplot2::geom_hline(yintercept = overall_sd, linetype = "dashed", color = "darkgreen") +
    ggplot2::labs(
      title    = "Binned Posterior Residual SD across Fitted Quantiles",
      subtitle = "Points = Posterior Median SD; Bars = 95% CI; Dashed = Overall SD",
      x        = "Fitted Value Quantile Bin (Low -> High)",
      y        = "Residual Standard Deviation"
    ) +
    ggplot2::theme_minimal()
  
  out <- list(
    binned_summary     = bin_summary,
    variance_ratio     = list(median = var_ratio_median, ci95 = var_ratio_ci),
    overall_sd         = overall_sd,
    plot_residuals     = p1,
    plot_binned        = p2,
    response           = resp
  )
  class(out) <- c("because_homoscedasticity", "list")
  return(out)
}

#' @export
print.because_homoscedasticity <- function(x, ...) {
  cat("\n=== Bayesian Homoscedasticity Diagnostics (Gelman & Hill 2007) ===\n")
  cat(paste("Response Variable:", x$response, "\n\n"))
  cat("Binned Residual Standard Deviations:\n")
  print(format(x$binned_summary[, c("Bin", "Median_SD", "Lower95", "Upper95")], justify = "left"), row.names = FALSE)
  
  cat(paste0("\nOverall Residual SD: ", round(x$overall_sd, 3), "\n"))
  cat(paste0("Variance Ratio (Top vs. Bottom Bin) Posterior Median: ", round(x$variance_ratio$median, 2), "\n"))
  cat(paste0("95% Credible Interval: [", round(x$variance_ratio$ci95[1], 2), ", ", round(x$variance_ratio$ci95[2], 2), "]\n"))
  cat("Interpretation: Credible intervals overlapping the overall SD line indicate homoscedasticity.\n\n")
  invisible(x)
}


# ==============================================================================
# 4. Residual Autocorrelation & Phylogenetic/Spatial Structure
# ==============================================================================

#' @title Bayesian Residual Structure & Variance Partition Diagnostics
#' @rdname check_residual_structure
#' @description
#' Evaluates the proportion of variance absorbed by phylogenetic or spatial random effects
#' in a \code{because} model (Hadfield 2010; de Villemereuil & Nakagawa 2014).
#'
#' Computes the posterior distribution of the **Variance Partition Coefficient (VPC)**:
#' \deqn{\lambda_{\text{Bayes}}^{(s)} = \frac{\sigma_{\text{struct}}^{2(s)}}{\sigma_{\text{struct}}^{2(s)} + \sigma_{\text{res}}^{2(s)}}}
#'
#' @param object A \code{because} model fit object.
#' @param resp Character string; name of the response variable. If \code{NULL}, takes the first response.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"because_residual_structure"}.
#'
#' @export
check_residual_structure.because <- function(object, resp = NULL, ...) {
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  samples_mat <- as.matrix(object$samples)
  
  # Search for structural standard deviations
  # Patterns: sigma_phylo_VAR, sigma_spatial_VAR, sigma_VAR_phylo, etc.
  struct_patterns <- paste0("sigma_(phylo|spatial|group)_", resp)
  struct_cols <- grep(struct_patterns, colnames(samples_mat), value = TRUE)
  if (length(struct_cols) == 0) {
    struct_cols <- grep(paste0("sigma_", resp, "_(phylo|spatial|group)"), colnames(samples_mat), value = TRUE)
  }
  
  if (length(struct_cols) == 0) {
    cat(paste("No phylogenetic or spatial covariance structure found for response:", resp, "\n"))
    return(invisible(NULL))
  }
  
  # Residual variance
  res_col <- paste0("sigma_", resp, "_res")
  if (!res_col %in% colnames(samples_mat)) {
    tau_col <- paste0("tau_res_", resp)
    if (tau_col %in% colnames(samples_mat)) {
      sigma_res <- 1 / sqrt(samples_mat[, tau_col])
    } else {
      sigma_res <- rep(1, nrow(samples_mat))
    }
  } else {
    sigma_res <- samples_mat[, res_col]
  }
  
  vpc_list <- list()
  for (sc in struct_cols) {
    sigma_struct <- samples_mat[, sc]
    vpc_samples <- (sigma_struct^2) / (sigma_struct^2 + sigma_res^2)
    struct_name <- gsub(paste0("sigma_|_", resp), "", sc)
    
    vpc_list[[struct_name]] <- data.frame(
      Structure = struct_name,
      Parameter = sc,
      Median    = round(stats::median(vpc_samples), 3),
      Lower95   = round(stats::quantile(vpc_samples, 0.025), 3),
      Upper95   = round(stats::quantile(vpc_samples, 0.975), 3),
      stringsAsFactors = FALSE
    )
  }
  
  vpc_summary <- do.call(rbind, vpc_list)
  
  out <- list(
    summary   = vpc_summary,
    response  = resp
  )
  class(out) <- c("because_residual_structure", "list")
  return(out)
}

#' @export
print.because_residual_structure <- function(x, ...) {
  cat("\n=== Bayesian Variance Partitioning / Heritability (Hadfield 2010) ===\n")
  cat(paste("Response Variable:", x$response, "\n\n"))
  cat("Posterior Distribution of Structure Variance Ratio (lambda_Bayes):\n")
  print(format(x$summary, justify = "left"), row.names = FALSE)
  cat("\nNote: lambda_Bayes = sigma_struct^2 / (sigma_struct^2 + sigma_res^2).\n")
  cat("Values approaching 1 indicate strong structural (phylogenetic/spatial) signal.\n\n")
  invisible(x)
}


# ==============================================================================
# 5. Unified Bayesian Diagnostic Dashboard
# ==============================================================================

#' @title Unified Bayesian Diagnostic Dashboard for because Models
#' @rdname check_model
#' @description
#' Generates a comprehensive multi-panel visual diagnostic summary for a \code{because} model equation.
#'
#' Panels included:
#' \enumerate{
#'   \item **Posterior Predictive Check** (\code{pp_check}): Observed density vs. replicated densities.
#'   \item **Bayesian Q-Q Plot** (\code{check_normality}): Randomized quantile residuals with credible simulation envelope.
#'   \item **Residuals vs. Fitted** (\code{check_homoscedasticity}): Pointwise residual spread.
#'   \item **Binned Residual Variance** (\code{check_homoscedasticity}): Variance dispersion across predictions.
#'   \item **Collinearity & Variance Inflation** (\code{check_collinearity}): Predictor BVIF and uncertainty inflation.
#' }
#'
#' @param object A \code{because} model fit object.
#' @param resp Character string; name of the response variable. If \code{NULL}, takes the first response.
#' @param ndraws Integer; number of draws to use for predictive checks. Defaults to 250.
#' @param ... Additional arguments.
#'
#' @return A list containing all diagnostic sub-objects and plots, printed invisibly.
#'
#' @importFrom graphics par layout
#' @export
check_model.because <- function(object, resp = NULL, ndraws = 250, ...) {
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  cat(paste("\nComputing Purely Bayesian Diagnostics for response:", resp, "...\n"))
  
  # 1. PPC
  p1 <- tryCatch(
    pp_check(object, resp = resp, ndraws = 50),
    error = function(e) NULL
  )
  
  # 2. Normality / Q-Q Envelope
  norm_res <- tryCatch(
    check_normality(object, resp = resp, ndraws = ndraws),
    error = function(e) NULL
  )
  p2 <- if (!is.null(norm_res)) norm_res$plot else NULL
  
  # 3. Homoscedasticity
  homo_res <- tryCatch(
    check_homoscedasticity(object, resp = resp, ndraws = ndraws),
    error = function(e) NULL
  )
  p3 <- if (!is.null(homo_res)) homo_res$plot_residuals else NULL
  p4 <- if (!is.null(homo_res)) homo_res$plot_binned else NULL
  
  # 4. Collinearity
  collin_res <- tryCatch(
    check_collinearity(object, equation = resp),
    error = function(e) NULL
  )
  p5 <- if (!is.null(collin_res)) plot(collin_res) else NULL
  
  # Assemble and display
  plots <- list(p1, p2, p3, p4, p5)
  valid_plots <- plots[!sapply(plots, is.null)]
  
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    do.call(gridExtra::grid.arrange, c(valid_plots, list(ncol = 2)))
  } else if (requireNamespace("patchwork", quietly = TRUE)) {
    # Combine with patchwork if installed
    combined <- patchwork::wrap_plots(valid_plots, ncol = 2)
    print(combined)
  } else {
    # Base grid layout fallback (always available in base R)
    n_p <- length(valid_plots)
    n_cols <- min(n_p, 2)
    n_rows <- ceiling(n_p / n_cols)
    
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(n_rows, n_cols)))
    
    for (i in seq_along(valid_plots)) {
      row_idx <- ceiling(i / n_cols)
      col_idx <- ((i - 1) %% n_cols) + 1
      print(valid_plots[[i]], vp = grid::viewport(layout.pos.row = row_idx, layout.pos.col = col_idx))
    }
  }
  
  out <- list(
    response         = resp,
    ppc_plot         = p1,
    normality        = norm_res,
    homoscedasticity = homo_res,
    collinearity     = collin_res
  )
  class(out) <- c("because_check_model", "list")
  invisible(out)
}
