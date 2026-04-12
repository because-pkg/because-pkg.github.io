#' Generate Posterior Predictive Draws for Because Models
#'
#' This function generates simulated data (\eqn{y_{rep}}) from the posterior distribution
#' of a fitted \code{because} model. These draws are used for posterior predictive checks (PPC)
#' and model validation.
#'
#' @param object A \code{because} fit object.
#' @param resp Character string; the name of the response variable to predict.
#'   If \code{NULL}, takes the first response variable in the model.
#' @param ndraws Integer; number of posterior draws to use. Defaults to all draws.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A matrix of dimensions \code{[ndraws x N_obs]} containing simulated response values.
#' @export
posterior_predict.because <- function(object, resp = NULL, ndraws = NULL, ...) {
  # 1. Selection & Identification
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  # 2. Extract Posterior Samples
  # Use as.matrix to combine chains
  samples_mat <- as.matrix(object$samples)
  n_total <- nrow(samples_mat)
  
  if (!is.null(ndraws) && ndraws < n_total) {
    # Systematic thinning to get ndraws
    set.seed(42)
    idx <- round(seq(1, n_total, length.out = ndraws))
    samples_mat <- samples_mat[idx, , drop = FALSE]
  }
  n_s <- nrow(samples_mat)
  
  # 3. Get Model Metadata
  pm <- object$parameter_map
  pm_resp <- pm[pm$response == resp, ]
  if (nrow(pm_resp) == 0) {
    stop(paste("Response variable", resp, "not found in model parameters."))
  }
  
  # 4. Reconstruct Linear Predictor (eta)
  # formula: eta = alpha + sum(beta * X) + sum(u_group[idx])
  
  # Get data resolution (N)
  data_list <- object$data
  # We use the original data length if it's a simple model, 
  # but for hierarchical we need the length of the response vector.
  obs_y <- object$original_data[[resp]]
  if (is.null(obs_y)) {
     # Try fetching from processed data list
     obs_y <- data_list[[resp]]
  }
  n_obs <- length(obs_y)
  
  # Use predict_ey helper logic (from marginal_effects) to build fixed part
  # We'll implement a simplified version here specialized for pp_check
  
  eta <- matrix(0, nrow = n_s, ncol = n_obs)
  
  # A. Intercept
  alpha_param <- pm_resp$parameter[pm_resp$predictor == "(Intercept)"]
  if (length(alpha_param) > 0 && alpha_param %in% colnames(samples_mat)) {
    eta <- eta + replicate(n_obs, samples_mat[, alpha_param])
  }
  
  # B. Fixed Effects (Slopes)
  beta_rows <- pm_resp[!pm_resp$predictor %in% c("(Intercept)", "(1 | Site)", "(1 | Survey)", "(1 | Species)"), ] # Basic check
  # Actually, need to exclude any random effect predictors
  beta_rows <- beta_rows[!grepl("\\|", beta_rows$predictor), ]
  
  for (i in seq_len(nrow(beta_rows))) {
    p_name <- beta_rows$parameter[i]
    v_name <- beta_rows$predictor[i]
    
    if (p_name %in% colnames(samples_mat) && v_name %in% names(object$original_data)) {
      b_samples <- samples_mat[, p_name]
      x_vals <- object$original_data[[v_name]]
      # Handle if x_vals is shorter (higher level predictor) by repeating it if id mapping is known
      # But for original data prediction, because usually has these as flat vectors already in JAGS data
      if (length(x_vals) != n_obs && v_name %in% names(data_list)) {
         x_vals <- data_list[[v_name]]
      }
      
      if (length(x_vals) == n_obs) {
        eta <- eta + (b_samples %*% t(as.matrix(x_vals)))
      }
    }
  }
  
  # C. Random Effects
  # Identify random effect nodes: u_resp_level[index]
  # We look for grep("^u_resp_", colnames(samples_mat))
  u_prefix <- paste0("u_", resp, "_")
  u_cols <- grep(paste0("^", u_prefix), colnames(samples_mat), value = TRUE)
  
  if (length(u_cols) > 0) {
    # Group by level name
    u_base_names <- unique(gsub("\\[\\d+\\]$", "", u_cols))
    
    for (ub in u_base_names) {
       lvl_name <- sub(u_prefix, "", ub)
       # Find the index mapping for this level in data_list
       # Usually matches "LevelName_idx"
       idx_name <- paste0(lvl_name, "_idx")
       if (!idx_name %in% names(data_list)) {
         # Search for any vector of length n_obs that looks like an index
         idx_name <- names(data_list)[sapply(data_list, function(x) length(x) == n_obs && max(x, na.rm=T) > 1)]
         # This is heuristic. Let's try harder to guess from the level mapping.
       }
       
       if (idx_name %in% names(data_list)) {
         indices <- data_list[[idx_name]]
         # Extract level samples into [n_s, J]
         sub_u_cols <- grep(paste0("^", ub, "\\["), colnames(samples_mat), value = TRUE)
         # Sort correctly
         sub_u_cols <- sub_u_cols[order(as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", sub_u_cols)))]
         u_samples <- samples_mat[, sub_u_cols, drop = FALSE]
         
         # Add for each observation
         eta <- eta + u_samples[, indices, drop = FALSE]
       }
    }
  }
  
  # 5. Inverse Link & Stochastic Draw
  fam <- object$family[[resp]] %||% "gaussian"
  
  y_rep <- matrix(0, nrow = n_s, ncol = n_obs)
  
  for (s in 1:n_s) {
    mu_s <- eta[s, ]
    
    if (fam == "gaussian") {
      # Need residual sigma
      sigma_name <- paste0("sigma_e_", resp)
      if (!sigma_name %in% colnames(samples_mat)) {
        tau_name <- paste0("tau_e_", resp)
        if (tau_name %in% colnames(samples_mat)) {
           sigma_val <- 1/sqrt(samples_mat[s, tau_name])
        } else {
           sigma_val <- 1 # Fixed
        }
      } else {
        sigma_val <- samples_mat[s, sigma_name]
      }
      y_rep[s, ] <- stats::rnorm(n_obs, mu_s, sigma_val)
      
    } else if (fam == "poisson") {
      y_rep[s, ] <- stats::rpois(n_obs, exp(mu_s))
      
    } else if (fam == "binomial") {
      prob_s <- 1 / (1 + exp(-mu_s))
      y_rep[s, ] <- stats::rbinom(n_obs, 1, prob_s)
    }
  }
  
  return(y_rep)
}

#' @export
posterior_predict <- function(object, ...) {
  UseMethod("posterior_predict")
}
