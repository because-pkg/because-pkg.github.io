#' @rdname posterior_predict
#'
#' @description
#' This function generates simulated data (\eqn{y_{rep}}) from the posterior distribution
#' of a fitted \code{because} model. These draws are used for posterior predictive checks (PPC)
#' and model validation.
#'
#' @param object A \code{because} fit object.
#' @param resp Character string; the name of the response variable to predict.
#'   If \code{NULL}, takes the first response variable in the model.
#' @param newdata Optional new data frame for prediction. Currently ignored.
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
      
      # Search top-level or sub-lists for predictor
      x_vals <- NULL
      search_data <- if (!is.null(newdata)) newdata else object$original_data
      
      if (v_name %in% names(search_data)) {
          x_vals <- search_data[[v_name]]
      } else if (is.list(search_data) && !is.data.frame(search_data)) {
          for (df_name in names(search_data)) {
              if (is.data.frame(search_data[[df_name]]) && v_name %in% names(search_data[[df_name]])) {
                  x_vals <- search_data[[df_name]][[v_name]]
                  break
              }
          }
      }
      
      if (is.null(x_vals)) x_vals <- data_list[[v_name]]

      if (is.null(x_vals)) {
          sanitized_v <- gsub(":", "_x_", v_name)
          if (sanitized_v %in% names(data_list)) {
              x_vals <- data_list[[sanitized_v]]
          } else if (sanitized_v %in% names(search_data)) {
              x_vals <- search_data[[sanitized_v]]
          } else if (grepl("(_x_|:)", v_name)) {
              parts <- strsplit(v_name, "(_x_|:)")[[1]]
              if (all(parts %in% names(search_data))) {
                  x_vals <- search_data[[parts[1]]] * search_data[[parts[2]]]
              } else if (all(parts %in% names(data_list))) {
                  x_vals <- data_list[[parts[1]]] * data_list[[parts[2]]]
              }
          }
      }
      
      # Hierarchical resolution alignment
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
              # Look for an index variable that maps these entities
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
  
  # Identify dispersion/variance parameter for Gaussian
  sigma_val <- rep(1, n_s)
  sigma_name <- paste0("sigma_", resp, "_res")
  tau_name <- paste0("tau_res_", resp)
  
  if (sigma_name %in% colnames(samples_mat)) {
      sigma_val <- samples_mat[, sigma_name]
  } else if (tau_name %in% colnames(samples_mat)) {
      sigma_val <- 1/sqrt(samples_mat[, tau_name])
  } else {
      # Fallback for unified or custom sigma names
      sig_match <- grep(paste0("^sigma_.*", resp, ".*res"), colnames(samples_mat), value = TRUE)
      if (length(sig_match) > 0) sigma_val <- samples_mat[, sig_match[1]]
  }

  # Identify dispersion parameter for Negative Binomial (r_<resp>)
  r_cols <- grep(paste0("^r_", resp, "($|_)"), colnames(samples_mat), value = TRUE)
  if (length(r_cols) == 0) {
      r_cols <- grep(paste0("^size_", resp, "($|_)"), colnames(samples_mat), value = TRUE)
  }

  # Identify zero-inflation probability for ZIP / ZINB (psi_<resp>)
  psi_cols <- grep(paste0("^psi_", resp, "($|_)"), colnames(samples_mat), value = TRUE)

  # Identify cutpoints for Ordinal
  cp_cols <- grep(paste0("^cutpoint_", resp, "(\\[|_|$)"), colnames(samples_mat), value = TRUE)
  if (length(cp_cols) > 0) {
      # Sort cutpoints in ascending order (cutpoint_Y[1], cutpoint_Y[2], ...)
      cp_nums <- as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", cp_cols))
      if (!any(is.na(cp_nums))) {
          cp_cols <- cp_cols[order(cp_nums)]
      }
  }

  for (s in 1:n_s) {
    mu_s <- eta[s, ]
    
    if (fam == "gaussian") {
      # If marginal, add the omitted RE variances to the simulation noise
      total_sd <- sqrt(sigma_val[s]^2 + re_variances[s])
      y_rep[s, ] <- stats::rnorm(n_obs, mean = mu_s, sd = total_sd)
    } else if (fam == "poisson") {
      clamped_mu <- pmin(20, pmax(-20, mu_s))
      y_rep[s, ] <- stats::rpois(n_obs, lambda = exp(clamped_mu))
    } else if (fam == "binomial") {
      clamped_mu <- pmin(20, pmax(-20, mu_s))
      prob_s <- stats::plogis(clamped_mu)
      y_rep[s, ] <- stats::rbinom(n_obs, size = 1, prob = prob_s)
    } else if (fam == "negbinomial") {
      clamped_mu <- pmin(20, pmax(-20, mu_s))
      r_val <- if (length(r_cols) > 0) samples_mat[s, r_cols[1]] else 1
      y_rep[s, ] <- stats::rnbinom(n_obs, size = r_val, mu = exp(clamped_mu))
    } else if (fam == "zip") {
      clamped_mu <- pmin(20, pmax(-20, mu_s))
      psi_val <- if (length(psi_cols) > 0) samples_mat[s, psi_cols[1]] else 0
      is_zero <- stats::rbinom(n_obs, size = 1, prob = psi_val)
      pois_draws <- stats::rpois(n_obs, lambda = exp(clamped_mu))
      y_rep[s, ] <- ifelse(is_zero == 1, 0, pois_draws)
    } else if (fam == "zinb") {
      clamped_mu <- pmin(20, pmax(-20, mu_s))
      psi_val <- if (length(psi_cols) > 0) samples_mat[s, psi_cols[1]] else 0
      r_val <- if (length(r_cols) > 0) samples_mat[s, r_cols[1]] else 1
      is_zero <- stats::rbinom(n_obs, size = 1, prob = psi_val)
      nb_draws <- stats::rnbinom(n_obs, size = r_val, mu = exp(clamped_mu))
      y_rep[s, ] <- ifelse(is_zero == 1, 0, nb_draws)
    } else if (fam == "ordinal" && length(cp_cols) > 0) {
      cps <- samples_mat[s, cp_cols]
      K <- length(cps) + 1
      for (i in 1:n_obs) {
        q_vals <- stats::plogis(cps - mu_s[i])
        probs <- numeric(K)
        probs[1] <- q_vals[1]
        if (K > 2) {
          for (k in 2:(K - 1)) {
            probs[k] <- max(0, q_vals[k] - q_vals[k - 1])
          }
        }
        probs[K] <- max(0, 1 - q_vals[K - 1])
        probs <- probs / sum(probs)
        y_rep[s, i] <- which(stats::rmultinom(1, size = 1, prob = probs) == 1)
      }
    } else {
      warning(sprintf("Posterior predictive simulation for family '%s' is not implemented. Returning NAs.", fam))
      y_rep[s, ] <- NA_real_
    }
  }
  
  return(y_rep)
}

#' Predict Method for because Model Fits
#'
#' Computes posterior predictions, fitted values, and credible intervals
#' for a fitted \code{because} model.
#'
#' @param object A \code{because} fit object.
#' @param resp Character string; name of the response variable to predict.
#'   If \code{NULL}, defaults to the first response variable.
#' @param newdata Optional new data frame for prediction. Currently ignored.
#' @param summary Logical; if \code{TRUE} (default), returns a data frame
#'   with summary statistics (\code{Estimate}, \code{Est.Error}, \code{Q2.5}, \code{Q97.5}).
#'   If \code{FALSE}, returns the raw draws matrix from \code{\link{posterior_predict}}.
#' @param probs Numeric vector of length 2; quantiles for summary intervals. Defaults to \code{c(0.025, 0.975)}.
#' @param ndraws Integer; number of posterior draws to use. Defaults to all draws.
#' @param re_formula Formula or \code{NA}; determines which random effects to include.
#'   See \code{\link{posterior_predict}} for details.
#' @param ... Additional arguments passed to \code{\link{posterior_predict}}.
#'
#' @return If \code{summary = TRUE}, a data frame containing:
#'   \item{Estimate}{Posterior mean prediction}
#'   \item{Est.Error}{Posterior standard deviation}
#'   \item{Q2.5}{Lower credible limit}
#'   \item{Q97.5}{Upper credible limit}
#'   If \code{summary = FALSE}, an \code{[ndraws x N_obs]} matrix of posterior draws.
#'
#' @export
predict.because <- function(
    object,
    resp = NULL,
    newdata = NULL,
    summary = TRUE,
    probs = c(0.025, 0.975),
    ndraws = NULL,
    re_formula = NULL,
    ...
) {
    draws <- posterior_predict(
        object,
        resp = resp,
        newdata = newdata,
        ndraws = ndraws,
        re_formula = re_formula,
        ...
    )
    if (!summary) {
        return(draws)
    }
    est <- colMeans(draws, na.rm = TRUE)
    se <- apply(draws, 2, stats::sd, na.rm = TRUE)
    q_lower <- apply(draws, 2, stats::quantile, probs = probs[1], na.rm = TRUE)
    q_upper <- apply(draws, 2, stats::quantile, probs = probs[2], na.rm = TRUE)

    res <- data.frame(
        Estimate = est,
        Est.Error = se
    )
    res[[paste0("Q", probs[1] * 100)]] <- q_lower
    res[[paste0("Q", probs[2] * 100)]] <- q_upper
    return(res)
}

# Helper %||% if not defined elsewhere
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @export
posterior_predict <- function(object, ...) {
  UseMethod("posterior_predict")
}
