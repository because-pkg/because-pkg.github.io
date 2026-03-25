#' Calculate Marginal Effects for Because Fit
#'
#' Estimates the Average Marginal Effect (AME) or Marginal Effect at the Mean (MEM)
#' for all paths in a model. This is especially useful for comparing coefficients across
#' different families (e.g., Gaussian vs. Binomial) as it converts them to a common scale
#' (expected change in the response).
#'
#' @param fit A because fit object.
#' @param at Character or list. If NULL (default), calculates Average Marginal Effects (AME).
#'   If "mean", calculates Marginal Effects at the Mean (MEM).
#' @param prob Numeric; probability mass for the credible interval (default 0.95).
#' @param samples Integer. Number of posterior samples to use (default 100 for speed).
#' @return A data frame with marginal effects per path.
#' @export
marginal_effects <- function(fit, at = NULL, prob = 0.95, samples = 100) {
  if (is.null(fit$parameter_map)) {
    stop("Fit object does not contain a parameter_map. Please refit the model.")
  }

  # 1. Extract posterior samples
  samples_mat <- as.matrix(fit$samples)
  n_total_samples <- nrow(samples_mat)
  
  # Subsample for performance
  if (!is.null(samples) && samples < n_total_samples) {
    idx <- round(seq(1, n_total_samples, length.out = samples))
    samples_mat <- samples_mat[idx, , drop = FALSE]
  }
  n_samples <- nrow(samples_mat)

  # 2. Identify all Response ~ Predictor paths
  # Exclude internal intercept/cutpoint "predictors"
  pm <- fit$parameter_map
  raw_paths <- pm[!pm$predictor %in% c("(Intercept)", "(Intercepts)", "(Cutpoints)"), ]
  
  # Group by response to avoid redundant calculations
  responses <- unique(raw_paths$response)
  
  results_list <- list()

  # 3. Process each response
  for (resp in responses) {
    resp_paths <- raw_paths[raw_paths$response == resp, ]
    
    # Safe family lookup
    dist <- "gaussian"
    if (!is.null(fit$family) && resp %in% names(fit$family)) {
       dist <- fit$family[[resp]]
    }
    
    # Get all predictors for this response (including intercepts)
    all_pm_resp <- pm[pm$response == resp, ]
    
    # Identify focal predictors (original variables)
    # We need to handle cases where a variable is expanded (Diet_L, Diet_Q)
    # We want the effect of "Diet".
    focal_vars <- unique(sapply(resp_paths$predictor, function(p) {
        # Check if it is a dummy/contract: var_L, var_Q, var_dummy, or var_k
        base <- gsub("(_L|_Q|_C|_dummy|_\\d+|\\[\\d+\\])$", "", p)
        return(base)
    }))

    # Get data for this response's level
    data_resp <- fit$original_data
    # (Handling hierarchical data would go here in future versions)
    
    # Define "baseline" individual(s)
    if (identical(at, "mean")) {
       # MEM: Marginal Effect at the Mean
       # Compute means of all predictors in pm
       base_data <- as.data.frame(lapply(data_resp, function(x) {
           if (is.numeric(x)) mean(x, na.rm=TRUE) else x[1]
       }))[1, , drop=FALSE]
    } else {
       # AME: Average Marginal Effect (over all observations)
       base_data <- data_resp
    }

    # Helper function to predict E[Y]
    # samples_chunk: matrix of samples [n_samples, p]
    # current_data: data.frame of predictors [N_obs, p]
    predict_ey <- function(samples_chunk, current_data, dist, response_name) {
       N_obs <- nrow(current_data)
       N_s <- nrow(samples_chunk)
       
       # Extract relevant parameters
       alpha_p <- all_pm_resp$parameter[all_pm_resp$predictor %in% c("(Intercept)", "(Intercepts)")]
       beta_rows <- all_pm_resp[!all_pm_resp$predictor %in% c("(Intercept)", "(Intercepts)", "(Cutpoints)"), ]
       cut_p <- all_pm_resp$parameter[all_pm_resp$predictor == "(Cutpoints)"]

       # Compute Linear Predictor Matrix [N_s, N_obs]
       # linpred = alpha + sum(beta_j * X_j)
       
       # Start with intercepts
       lp_mat <- matrix(0, nrow = N_s, ncol = N_obs)
       
       # Add Alpha
       if (length(alpha_p) > 0 && alpha_p %in% colnames(samples_chunk)) {
          # alpha_p might be a vector for multinomial? 
          # No, for multinomial we handle it differently.
          lp_mat <- lp_mat + replicate(N_obs, samples_chunk[, alpha_p])
       }

       # Add Betas
       for(i in seq_len(nrow(beta_rows))) {
          p_name <- beta_rows$parameter[i]
          v_name <- beta_rows$predictor[i]
          
          if (p_name %in% colnames(samples_chunk) && v_name %in% colnames(current_data)) {
             beta_s <- samples_chunk[, p_name] # [N_s]
             x_v <- current_data[[v_name]]     # [N_obs]
             lp_mat <- lp_mat + (beta_s %*% t(as.matrix(x_v)))
          }
       }
       
       # Convert Linear Predictor to E[Y]
       if (dist == "gaussian") {
          return(lp_mat)
       } else if (dist == "binomial") {
          return(1 / (1 + exp(-lp_mat))) # Logit -> Prob
       } else if (dist == "ordinal") {
          # Ordinal Expected Value: sum(k * p_k)
          # We need cutpoints
          # Param map store cutpoint as "cutpoint_var"
          # samples has "cutpoint_var[1]", "cutpoint_var[2]" ...
          base_cut_name <- gsub("\\[\\]$", "", cut_p)
          all_param_names <- colnames(samples_chunk)
          cut_params <- grep(paste0("^", base_cut_name, "\\[\\d+\\]$"), all_param_names, value=TRUE)
          # Sort them 1, 2, 3
          cut_params <- cut_params[order(as.numeric(gsub(".*\\[(\\d+)\\]$", "\\1", cut_params)))]
          
          K <- length(cut_params) + 1
          probs_list <- list()
          
          # Cumulative probabilities P(Y <= k) = logit_inv(cut_k - lp)
          prev_phi <- matrix(0, nrow=N_s, ncol=N_obs)
          ey_mat <- matrix(0, nrow=N_s, ncol=N_obs)
          
          for (k in 1:(K-1)) {
             cut_s <- samples_chunk[, cut_params[k]] # [N_s]
             # Phi_k = logit_inv(cut_k - lp)
             phi_k <- 1 / (1 + exp(-(replicate(N_obs, cut_s) - lp_mat)))
             p_k <- phi_k - prev_phi
             ey_mat <- ey_mat + k * p_k
             prev_phi <- phi_k
          }
          # Last category K
          p_K <- 1 - prev_phi
          ey_mat <- ey_mat + K * p_K
          return(ey_mat)
       } else {
          # Fallback or unknown
          return(lp_mat)
       }
    }

    # 4. Calculate Marginal Effect for each focal variable
    for (f_var in focal_vars) {
       # Base Prediction [N_s, N_obs]
       ey_base <- predict_ey(samples_mat, base_data, dist, resp)
       
       # Plus One Prediction
       plus_data <- base_data
       if (f_var %in% names(plus_data)) {
          # Increment Focal Variable
          if (is.numeric(plus_data[[f_var]])) {
             plus_data[[f_var]] <- plus_data[[f_var]] + 1
          }
          
          # Update Derived Polynomial Terms
          if (!is.null(fit$poly_terms)) {
             for (pt in fit$poly_terms) {
                if (pt$base_var == f_var) {
                   plus_data[[pt$internal_name]] <- plus_data[[f_var]]^pt$power
                }
             }
          }
          
          # Update Categorical Derived Terms (Dummies/Contrasts)
          if (!is.null(fit$categorical_vars) && f_var %in% names(fit$categorical_vars)) {
             # This is tricky without the full contrast matrix.
             # In v1, we assume numeric incrementing of the contrasts themselves,
             # which is an approximation for 'a unit change in the factor latent scale'.
             # More robust handling would require re-running contr.poly on the new levels.
          }
       }
       
       ey_plus <- predict_ey(samples_mat, plus_data, dist, resp)
       
       # Difference
       me_samples <- rowMeans(ey_plus - ey_base, na.rm=TRUE) # [N_s]
       
       # Summarize
       results_list[[length(results_list) + 1]] <- data.frame(
          Response = resp,
          Predictor = f_var,
          Effect = mean(me_samples, na.rm=TRUE),
          Lower = quantile(me_samples, (1 - prob)/2, na.rm=TRUE),
          Upper = quantile(me_samples, 1 - (1 - prob)/2, na.rm=TRUE),
          Family = dist,
          stringsAsFactors = FALSE
       )
    }
  }

  res_df <- do.call(rbind, results_list)
  rownames(res_df) <- NULL
  return(res_df)
}
