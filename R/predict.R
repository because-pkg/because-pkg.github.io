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
#' @param re_formula Formula or \code{NA}; determines which random effects to include.
#'   If \code{NULL} (default), all random effects are included (conditional prediction).
#'   If \code{NA}, no random effects are included (marginal prediction).
#' @param ... Additional arguments (currently ignored).
#'
#' @return A matrix of dimensions \code{[ndraws x N_obs]} containing simulated response values.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(Y = rnorm(100), X = rnorm(100))
#' fit <- because(list(Y ~ X), data = df)
#'
#' # Draw 200 posterior predictive samples for Y
#' yrep <- posterior_predict(fit, resp = "Y", ndraws = 200)
#' dim(yrep)  # [200 x 100]
#'
#' # Marginal prediction (no random effects)
#' yrep_marg <- posterior_predict(fit, resp = "Y", ndraws = 200, re_formula = NA)
#' }
#' @export
posterior_predict.because <- function(object, resp = NULL, newdata = NULL, ndraws = NULL, re_formula = NULL, ...) {
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
  obs_y <- if (!is.null(newdata)) {
     if (is.data.frame(newdata)) newdata[[resp]] else if (is.list(newdata) && !is.data.frame(newdata)) {
         # Find the data frame containing the response
         found <- NULL
         for (df_name in names(newdata)) {
             if (is.data.frame(newdata[[df_name]]) && resp %in% names(newdata[[df_name]])) {
                 found <- newdata[[df_name]][[resp]]
                 break
             }
         }
         found
     } else newdata[[resp]]
  } else {
     object$original_data[[resp]]
  }
  
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

  n_obs <- if (is.matrix(obs_y)) ncol(obs_y) else length(obs_y)
  if (target_level != "obs" && !is.null(h_info$link_vars[[target_level]])) {
      idx_var <- h_info$link_vars[[target_level]]
      # Search original_data or data for unique entities
      search_data <- if (!is.null(newdata)) newdata else object$original_data
      if (idx_var %in% names(search_data)) {
          n_obs <- length(unique(search_data[[idx_var]]))
      } else if (paste0(idx_var, "_idx") %in% names(data_list)) {
          n_obs <- max(data_list[[paste0(idx_var, "_idx")]], na.rm=TRUE)
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
  beta_rows <- pm_resp[pm_resp$type == "coefficient" & pm_resp$predictor != "(Intercept)", ]
  
  for (i in seq_len(nrow(beta_rows))) {
    p_name <- beta_rows$parameter[i]
    v_name <- beta_rows$predictor[i]
    
    if (p_name %in% colnames(samples_mat)) {
      b_samples <- samples_mat[, p_name]
      
      # [Smart Lookup] Search top-level or sub-lists for predictor
      x_vals <- NULL
      search_data <- if (!is.null(newdata)) newdata else object$original_data
      
      if (v_name %in% names(search_data)) {
          x_vals <- search_data[[v_name]]
      } else if (is.list(search_data) && !is.data.frame(search_data)) {
          # Search inside site, survey, etc.
          for (df_name in names(search_data)) {
              if (is.data.frame(search_data[[df_name]]) && v_name %in% names(search_data[[df_name]])) {
                  x_vals <- search_data[[df_name]][[v_name]]
                  break
              }
          }
      }
      
      # If still not found, check standard data list
      if (is.null(x_vals)) x_vals <- data_list[[v_name]]
      
      # [Hierarchical Resolution Alignment]
      if (!is.null(x_vals)) {
          # Case 1: x_vals is observation-level but we need higher level (e.g. species traits)
          # Check length of x_vals; if it's a matrix, check ncol
          len_x <- if (is.matrix(x_vals)) ncol(x_vals) else length(x_vals)
          
          if (len_x > n_obs && target_level != "obs") {
              idx_var <- h_info$link_vars[[target_level]]
              # Search for links in search_data
              links <- NULL
              for (df_name in names(search_data)) {
                  if (is.data.frame(search_data[[df_name]]) && idx_var %in% names(search_data[[df_name]])) {
                      links <- search_data[[df_name]][[idx_var]]
                      # Ensure links has same length as x_vals
                      if (length(links) == len_x) break else links <- NULL
                  }
              }
              if (!is.null(links)) {
                  unique_idx <- !duplicated(links)
                  if (is.matrix(x_vals)) x_vals <- x_vals[, unique_idx, drop = FALSE] else x_vals <- x_vals[unique_idx]
                  len_x <- sum(unique_idx)
              }
          }
          
          # Case 2: x_vals is site-level but response is survey-level (Elevation_s -> Temperature)
          # We need to map len_x to n_obs
          if (len_x < n_obs) {
              # Look for a JAGS index variable that maps these entities
              idx_candidate <- NULL
              for (idx_name in grep("_idx$", names(data_list), value = TRUE)) {
                  if (length(data_list[[idx_name]]) == n_obs && max(data_list[[idx_name]], na.rm = TRUE) == len_x) {
                      idx_candidate <- data_list[[idx_name]]
                      break
                  }
              }
              
              if (!is.null(idx_candidate)) {
                  if (is.matrix(x_vals)) {
                      x_vals <- x_vals[, idx_candidate, drop = FALSE]
                  } else {
                      x_vals <- x_vals[idx_candidate]
                  }
                  len_x <- n_obs
              } else if (n_obs %% len_x == 0) {
                  # Fallback for perfectly balanced implicit hierarchies
                  rep_times <- n_obs / len_x
                  if (is.matrix(x_vals)) {
                      x_vals <- x_vals[, rep(1:len_x, each = rep_times), drop = FALSE]
                  } else {
                      x_vals <- rep(x_vals, each = rep_times)
                  }
                  len_x <- n_obs
              }
          }
      }
      
      if (!is.null(x_vals) && len_x == n_obs) {
        contribution <- if (is.matrix(x_vals)) {
            (as.numeric(b_samples) * x_vals)
        } else {
            (b_samples %*% t(as.matrix(x_vals)))
        }
        eta <- eta + contribution
      }
    }
  }
  
  message("DEBUG [", resp, "]: eta range: ", paste(round(range(eta), 3), collapse = " to "))
  # We track RE variances for marginal checks (re_formula = NA)
  re_variances <- matrix(0, nrow = n_s, ncol = 1)

  u_prefix <- paste0("u_", resp, "_")
  # Determine all possible random effect variances for this response
  sigma_re_names <- grep(paste0("^sigma_.*_", resp, "$"), colnames(samples_mat), value = TRUE)
  # Also catch sigma_phylo_VAR, sigma_site_VAR etc.
  sigma_re_names <- unique(c(sigma_re_names, grep(paste0("^sigma_[^res].*_", resp, "($|_)"), colnames(samples_mat), value = TRUE)))

  if (is.null(re_formula) || !is.na(re_formula)) {
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
                   eta <- eta + u_samples[, 1:n_obs, drop = FALSE]
               }
            }
          }
      }
    }
  } else if (!is.null(re_formula) && is.na(re_formula)) {
      # Marginal Prediction: Collect RE variances to add to the predictive noise
      for (sn in sigma_re_names) {
          re_variances <- re_variances + samples_mat[, sn, drop=FALSE]^2
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
  }
  
  for (s in 1:n_s) {
    mu_s <- eta[s, ]
    
    if (fam == "gaussian") {
      # If marginal, add the omitted RE variances to the simulation noise
      total_sd <- sqrt(sigma_val[s]^2 + re_variances[s])
      y_rep[s, ] <- stats::rnorm(n_obs, mu_s, total_sd)
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
