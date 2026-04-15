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
  
  # 3. Get Model Metadata & Determine Resolution
  pm <- object$parameter_map
  pm_resp <- pm[pm$response == resp, ]
  if (nrow(pm_resp) == 0) {
    stop(paste("Response variable", resp, "not found in model parameters."))
  }
  
  data_list <- object$data
  obs_y <- object$original_data[[resp]]
  if (is.null(obs_y)) {
     obs_y <- data_list[[resp]]
  }
  
  # Identify hierarchical level of the response
  h_info <- object$hierarchical_info
  target_level <- "obs"
  if (!is.null(h_info)) {
    for (lvl in names(h_info$levels)) {
       if (resp %in% h_info$levels[[lvl]]) {
          target_level <- lvl
          break
       }
    }
  }

  n_obs <- length(obs_y)
  if (target_level != "obs" && !is.null(h_info$link_vars[[target_level]])) {
      idx_var <- h_info$link_vars[[target_level]]
      # Search original_data or data for unique entities
      if (idx_var %in% names(object$original_data)) {
          n_obs <- length(unique(object$original_data[[idx_var]]))
      } else if (paste0(idx_var, "_idx") %in% names(data_list)) {
          n_obs <- max(data_list[[paste0(idx_var, "_idx")]], na.rm=T)
      }
  }
  
  # 4. Reconstruct Linear Predictor (eta)
  eta <- matrix(0, nrow = n_s, ncol = n_obs)
  
  # A. Intercept
  alpha_param <- pm_resp$parameter[pm_resp$predictor == "(Intercept)"]
  if (length(alpha_param) > 0 && alpha_param %in% colnames(samples_mat)) {
    eta <- eta + replicate(n_obs, samples_mat[, alpha_param])
  }
  
  # B. Fixed Effects (Slopes)
  # Basic check excluding intercept and random effects
  beta_rows <- pm_resp[!pm_resp$predictor %in% c("(Intercept)"), ]
  beta_rows <- beta_rows[!grepl("\\|", beta_rows$predictor) & !grepl("^u_", beta_rows$parameter), ]
  
  for (i in seq_len(nrow(beta_rows))) {
    p_name <- beta_rows$parameter[i]
    v_name <- beta_rows$predictor[i]
    
    if (p_name %in% colnames(samples_mat)) {
      b_samples <- samples_mat[, p_name]
      # Predictor prioritizes the target level resolution
      x_vals <- object$original_data[[v_name]]
      
      # If x_vals is observation-level but we are predicting at higher level, subset it
      if (!is.null(x_vals) && length(x_vals) > n_obs && target_level != "obs") {
          idx_var <- h_info$link_vars[[target_level]]
          if (!is.null(idx_var) && idx_var %in% names(object$original_data)) {
              links <- object$original_data[[idx_var]]
              # Take values for the unique entities only
              x_vals <- x_vals[!duplicated(links)]
          }
      }
      
      # Fallback to data_list if missing/wrong length
      if (is.null(x_vals) || length(x_vals) != n_obs) {
         if (v_name %in% names(data_list)) {
             x_v <- data_list[[v_name]]
             if (length(x_v) == n_obs) x_vals <- x_v
         }
      }
      
      if (!is.null(x_vals) && length(x_vals) == n_obs) {
        eta <- eta + (b_samples %*% t(as.matrix(x_vals)))
      }
    }
  }
  
  # C. Random Effects & Structures
  u_prefix <- paste0("u_", resp, "_")
  u_cols <- grep(paste0("^", u_prefix), colnames(samples_mat), value = TRUE)
  
  if (length(u_cols) > 0) {
    u_base_names <- unique(gsub("\\[\\d+\\]$", "", u_cols))
    
    for (ub in u_base_names) {
        lvl_name <- sub(u_prefix, "", ub)
        
        # Mapping index
        if (lvl_name == target_level || (lvl_name == "phylo" && target_level == "species")) {
            indices <- 1:n_obs
        } else {
            idx_name <- NULL
            if (!is.null(h_info) && !is.null(h_info$link_vars)) {
                 if (lvl_name %in% names(h_info$link_vars)) {
                     idx_var_name <- h_info$link_vars[[lvl_name]]
                     idx_name <- paste0(idx_var_name, "_idx")
                 } else if (lvl_name == "phylo" && "species" %in% names(h_info$link_vars)) {
                     idx_var_name <- h_info$link_vars[["species"]]
                     idx_name <- paste0(idx_var_name, "_idx")
                 }
            }
            if (is.null(idx_name) || !idx_name %in% names(data_list)) {
                idx_name <- paste0(lvl_name, "_idx")
            }
            indices <- if (idx_name %in% names(data_list)) data_list[[idx_name]] else NULL
        }
        
        if (!is.null(indices)) {
          # Extract and sort level samples
          sub_u_cols <- grep(paste0("^", ub, "\\["), colnames(samples_mat), value = TRUE)
          indices_numeric <- as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", sub_u_cols))
          sub_u_cols <- sub_u_cols[order(indices_numeric)]
          u_samples <- samples_mat[, sub_u_cols, drop = FALSE]
          
          valid_idx <- !is.na(indices)
          safe_indices <- indices[valid_idx]
          
          if (length(safe_indices) > 0) {
             if (length(indices) == ncol(eta)) {
                eta[, valid_idx] <- eta[, valid_idx] + u_samples[, safe_indices, drop = FALSE]
             } else if (length(indices) > ncol(eta) && target_level != "obs") {
                 # Fallback: if we only have one-per-entity prediction but long indices, take first J
                 eta <- eta + u_samples[, 1:n_obs, drop = FALSE]
             }
          }
        }
    }
  }
  
  # 5. Inverse Link & Stochastic Draw
  fam <- if (!is.null(object$family) && resp %in% names(object$family)) object$family[[resp]] else "gaussian"
  y_rep <- matrix(0, nrow = n_s, ncol = n_obs)
  
  # Identify dispersion/variance parameter
  sigma_val <- rep(1, n_s)
  sigma_name <- paste0("sigma_", resp, "_res")
  tau_name <- paste0("tau_res_", resp)
  
  if (sigma_name %in% colnames(samples_mat)) {
      sigma_val <- samples_mat[, sigma_name]
  } else if (tau_name %in% colnames(samples_mat)) {
      sigma_val <- 1/sqrt(samples_mat[, tau_name])
  } else {
      # Fallback to legacy naming for backward compatibility
      sigma_legacy <- paste0("sigma_e_", resp)
      tau_legacy <- paste0("tau_e_", resp)
      if (sigma_legacy %in% colnames(samples_mat)) {
          sigma_val <- samples_mat[, sigma_legacy]
      } else if (tau_legacy %in% colnames(samples_mat)) {
          sigma_val <- 1/sqrt(samples_mat[, tau_legacy])
      }
  }
  
  for (s in 1:n_s) {
    mu_s <- eta[s, ]
    
    if (fam == "gaussian") {
      y_rep[s, ] <- stats::rnorm(n_obs, mu_s, sigma_val[s])
    } else if (fam == "poisson") {
      y_rep[s, ] <- stats::rpois(n_obs, exp(mu_s))
    } else if (fam == "binomial") {
      prob_s <- 1 / (1 + exp(-mu_s))
      y_rep[s, ] <- stats::rbinom(n_obs, 1, prob_s)
    } else if (fam == "negbinomial") {
      size_name <- paste0("size_", resp)
      cur_size <- if (size_name %in% colnames(samples_mat)) samples_mat[s, size_name] else 1
      y_rep[s, ] <- stats::rnbinom(n_obs, size = cur_size, mu = exp(mu_s))
    }
  }
  
  return(y_rep)
}

# Helper %||% if not defined elsewhere
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @export
posterior_predict <- function(object, ...) {
  UseMethod("posterior_predict")
}
