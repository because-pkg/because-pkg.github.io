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
#' @param multinomial_probabilities Logical. If TRUE, returns granular probability shifts for each category
#'   of multinomial (unordered) responses instead of a single expected value shift. Default FALSE.
#' @return A data frame with marginal effects per path.
#' @export
marginal_effects <- function(fit, at = NULL, prob = 0.95, samples = 100, multinomial_probabilities = FALSE) {
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

  # Helper for robust parameter matching (handles var[] and var[k])
  get_param_samples <- function(base_names, samples) {
    if (length(base_names) == 0) return(NULL)
    
    all_matches <- c()
    cols <- colnames(samples)
    
    for (bn in base_names) {
       clean_name <- gsub("\\[\\]$", "", bn)
       
       # Exact match
       if (clean_name %in% cols) {
          all_matches <- c(all_matches, clean_name)
       } else {
          # Array match
          matches <- grep(paste0("^", clean_name, "\\[\\d+\\]$"), cols, value=TRUE)
          if (length(matches) > 0) {
             # Sort matches to ensure [1], [2], [3] order
             matches <- matches[order(as.numeric(gsub(".*\\[(\\d+)\\].*", "\\1", matches)))]
             all_matches <- c(all_matches, matches)
          }
       }
    }
    
    if (length(all_matches) > 0) {
       # Return as matrix
       return(samples[, unique(all_matches), drop=FALSE])
    }
    return(NULL)
  }

  # Helper function to predict E[Y]
  # samples_chunk: matrix of samples [n_samples, p]
  # current_data: data.frame of predictors [N_obs, p]
  predict_ey <- function(samples_chunk, current_data, dist, response_name, all_pm_resp, multinomial_probabilities = FALSE) {
    N_obs <- nrow(current_data)
    N_s <- nrow(samples_chunk)
    
    if (is.null(N_obs) || is.null(N_s) || is.na(N_obs) || is.na(N_s) || !is.numeric(N_obs) || !is.numeric(N_s)) {
      stop(sprintf("Failed to determine dimensions for response '%s' (N_s=%s, N_obs=%s).", 
           response_name, as.character(N_s %||% "NULL"), as.character(N_obs %||% "NULL")))
    }

    # Extract relevant parameters
    alpha_p <- all_pm_resp$parameter[all_pm_resp$predictor %in% c("(Intercept)", "(Intercepts)")]
    beta_rows <- all_pm_resp[!all_pm_resp$predictor %in% c("(Intercept)", "(Intercepts)", "(Cutpoints)"), ]
    cut_p <- all_pm_resp$parameter[all_pm_resp$predictor == "(Cutpoints)"]

    # Compute Linear Predictor Matrix [N_s, N_obs]
    # linpred = alpha + sum(beta_j * X_j)
    
    # Start with intercepts
    lp_mat <- matrix(0, nrow = N_s, ncol = N_obs)
    
    if (dist == "multinomial") {
      # Multinomial Logic: Identify K and extract parameters per category
      # Get K from categorical_vars if available
      cat_info_resp <- fit$categorical_vars[[response_name]]
      K <- if (!is.null(cat_info_resp)) length(cat_info_resp$levels) else {
        # Fallback: Detect K from alpha or beta array dimensions
        alpha_p_array <- alpha_p[grepl("\\[\\]$", alpha_p)]
        temp_s <- if(length(alpha_p_array) > 0) get_param_samples(alpha_p_array, samples_chunk) else NULL
        if (!is.null(temp_s)) ncol(temp_s) + 1 else 2
      }
      
      # Initialize Probability Array [N_s, N_obs, K]
      # Level 1 is ALWAYS the reference (score = 0)
      scores <- array(0, dim=c(N_s, N_obs, K)) 
      
      base_alpha <- if(length(alpha_p) > 0) gsub("\\[\\]$", "", alpha_p[1]) else NULL
      
      for (k in 2:K) {
        # 1. Intercept for category k
        # We try alpha[k] (if alpha[1] exists) or alpha[k-1] (if ref is skipped)
        p_name <- paste0(base_alpha, "[", k, "]")
        if (p_name %in% colnames(samples_chunk)) {
          scores[,,k] <- scores[,,k] + replicate(N_obs, samples_chunk[, p_name])
        } else {
          # Try shift index k-1
          p_name_shift <- paste0(base_alpha, "[", k-1, "]")
          if (p_name_shift %in% colnames(samples_chunk)) {
            scores[,,k] <- scores[,,k] + replicate(N_obs, samples_chunk[, p_name_shift])
          }
        }
        
        # 2. Slopes for category k
        for(i in seq_len(nrow(beta_rows))) {
          v_name <- beta_rows$predictor[i]
          base_beta <- gsub("\\[\\]$", "", beta_rows$parameter[i])
          
          # Try beta[k] or beta[k-1]
          bp_name <- paste0(base_beta, "[", k, "]")
          if (!bp_name %in% colnames(samples_chunk)) {
             bp_name <- paste0(base_beta, "[", k-1, "]")
          }
          
          if (bp_name %in% colnames(samples_chunk) && v_name %in% colnames(current_data)) {
            beta_s <- samples_chunk[, bp_name]
            x_v <- current_data[[v_name]]
            scores[,,k] <- scores[,,k] + (beta_s %*% t(as.matrix(x_v)))
          }
        }
      }
      
      # Softmax to get probabilities
      exp_scores <- exp(scores)
      sum_exp <- apply(exp_scores, c(1,2), sum)
      
      if (multinomial_probabilities) {
        probs <- array(0, dim=c(N_s, N_obs, K))
        for (k in 1:K) {
          probs[,,k] <- exp_scores[,,k] / sum_exp
        }
        return(probs)
      } else {
        ey_mat <- matrix(0, nrow=N_s, ncol=N_obs)
        for (k in 1:K) {
          p_k <- exp_scores[,,k] / sum_exp
          ey_mat <- ey_mat + k * p_k
        }
        return(ey_mat)
      }
      
    } else {
      # Gaussian / Binomial / Ordinal Logic
      
      # Add Intercept
      if (length(alpha_p) > 0) {
        alpha_s <- get_param_samples(alpha_p, samples_chunk)
        if (!is.null(alpha_s)) {
          lp_mat <- lp_mat + replicate(N_obs, alpha_s[, 1])
        }
      }

      # Add Slopes
      for(i in seq_len(nrow(beta_rows))) {
        p_name <- beta_rows$parameter[i]
        v_name <- beta_rows$predictor[i]
        
        beta_s <- get_param_samples(p_name, samples_chunk)
        if (!is.null(beta_s) && v_name %in% colnames(current_data)) {
          lp_mat <- lp_mat + (beta_s[,1] %*% t(as.matrix(current_data[[v_name]])))
        }
      }
      
      if (dist == "gaussian") {
        return(lp_mat)
      } else if (dist == "binomial") {
        return(1 / (1 + exp(-lp_mat)))
      } else if (dist == "ordinal") {
        base_cut_name <- gsub("\\[\\]$", "", cut_p)
        cut_samples <- get_param_samples(base_cut_name, samples_chunk)
        if (is.null(cut_samples)) return(lp_mat)

        K <- ncol(cut_samples) + 1
        prev_phi <- matrix(0, nrow=N_s, ncol=N_obs)
        ey_mat <- matrix(0, nrow=N_s, ncol=N_obs)
        
        for (k in 1:(K-1)) {
          phi_k <- 1 / (1 + exp(-(replicate(N_obs, cut_samples[, k]) - lp_mat)))
          p_k <- phi_k - prev_phi
          ey_mat <- ey_mat + k * p_k
          prev_phi <- phi_k
        }
        ey_mat <- ey_mat + K * (1 - prev_phi)
        return(ey_mat)
      } else {
        return(lp_mat)
      }
    }
  }

  # 3. Process each response
  for (resp in responses) {
    resp_paths <- raw_paths[raw_paths$response == resp, ]
    
    dist <- "gaussian"
    if (!is.null(fit$family) && resp %in% names(fit$family)) {
       dist <- fit$family[[resp]]
    }
    
    all_pm_resp <- pm[pm$response == resp, ]
    
    focal_vars <- unique(sapply(resp_paths$predictor, function(p) {
        gsub("(_L|_Q|_C|_dummy|_\\d+|\\[\\d+\\])$", "", p)
    }))

    data_resp <- fit$original_data
    if (is.null(data_resp)) {
       if (is.data.frame(fit$data)) {
          data_resp <- fit$data
       } else if (is.list(fit$data)) {
          lengths <- sapply(fit$data, length)
          N_guess <- as.numeric(names(sort(table(lengths), decreasing=TRUE))[1])
          valid_cols <- names(lengths)[lengths == N_guess]
          data_resp <- as.data.frame(fit$data[valid_cols])
       }
    }
    
    if (is.null(data_resp) || nrow(data_resp) == 0) {
       stop("Missing data in fit object.")
    }

    if (identical(at, "mean")) {
       base_data <- as.data.frame(lapply(data_resp, function(x) {
           if (is.numeric(x)) mean(x, na.rm=TRUE) else x[1]
       }))[1, , drop=FALSE]
    } else {
       base_data <- data_resp
    }

    # 4. Calculate Marginal Effect for each focal variable
    for (f_var in focal_vars) {
       # Check if the predictor is an unordered multinomial (factor) variable
       pred_cat_info <- fit$categorical_vars[[f_var]]
       pred_is_multinomial <- !is.null(pred_cat_info) && pred_cat_info$type != "ordered"

       if (pred_is_multinomial && multinomial_probabilities) {
         # --- Multinomial PREDICTOR: compare each category k vs. reference (k=1) ---
         K_pred        <- length(pred_cat_info$levels)
         pred_levels   <- pred_cat_info$levels

         # Build the reference (category 1) data for all observations
         ref_data <- base_data
         for (i in 2:K_pred) {
           d_name <- paste0(f_var, "_", pred_cat_info$levels[i])
           if (d_name %in% names(ref_data)) ref_data[[d_name]] <- 0L
         }
         ey_ref <- predict_ey(samples_mat, ref_data, dist, resp, all_pm_resp, FALSE)

         for (k in 2:K_pred) {
           # Set data to category k
           cat_data <- ref_data
           d_name   <- paste0(f_var, "_", pred_cat_info$levels[k])
           if (d_name %in% names(cat_data)) cat_data[[d_name]] <- 1L

           ey_cat      <- predict_ey(samples_mat, cat_data, dist, resp, all_pm_resp, FALSE)
           me_samples  <- rowMeans(ey_cat - ey_ref, na.rm = TRUE)

           results_list[[length(results_list) + 1]] <- data.frame(
             Response         = resp,
             Predictor        = f_var,
             Category         = pred_levels[k],
             Effect           = mean(me_samples, na.rm = TRUE),
             Lower            = quantile(me_samples, (1 - prob) / 2, na.rm = TRUE),
             Upper            = quantile(me_samples, 1 - (1 - prob) / 2, na.rm = TRUE),
             Family           = dist,
             stringsAsFactors = FALSE
           )
         }

       } else {
         # --- Standard predictor: base prediction vs. +1 increment ---
         ey_base <- predict_ey(samples_mat, base_data, dist, resp, all_pm_resp, multinomial_probabilities)

         plus_data <- base_data
         if (f_var %in% names(plus_data)) {
           if (is.numeric(plus_data[[f_var]])) {
             plus_data[[f_var]] <- plus_data[[f_var]] + 1
           }

           if (!is.null(fit$poly_terms)) {
             for (pt in fit$poly_terms) {
               if (pt$base_var == f_var) {
                 plus_data[[pt$internal_name]] <- plus_data[[f_var]]^pt$power
               }
             }
           }

           if (!is.null(fit$categorical_vars) && f_var %in% names(fit$categorical_vars)) {
             cat_info <- fit$categorical_vars[[f_var]]
             K <- length(cat_info$levels)
             curr_idx <- as.integer(plus_data[[f_var]])

             if (cat_info$type == "ordered") {
               c_mat    <- stats::contr.poly(K)
               idx_plus <- pmin(curr_idx + 1, K)
               for (i in seq_along(cat_info$dummies)) {
                 d_name <- cat_info$dummies[i]
                 if (d_name %in% names(plus_data)) {
                   plus_data[[d_name]] <- c_mat[idx_plus, i]
                 }
               }
             } else {
               idx_plus <- pmin(curr_idx + 1, K)
               for (i in 2:K) {
                 d_name <- paste0(f_var, "_", cat_info$levels[i])
                 if (d_name %in% names(plus_data)) {
                   plus_data[[d_name]] <- as.numeric(idx_plus == i)
                 }
               }
             }
           }
         }

         ey_plus <- predict_ey(samples_mat, plus_data, dist, resp, all_pm_resp, multinomial_probabilities)

         # Aggregate Results
         if (dist == "multinomial" && multinomial_probabilities) {
           K             <- dim(ey_base)[3]
           resp_cat_info <- fit$categorical_vars[[resp]]
           resp_levels   <- if (!is.null(resp_cat_info)) resp_cat_info$levels else paste0("Level", 1:K)

           for (k in 1:6) { # Note: Hardcoded limits like this are risky, but K was 1:K previously. Fixed below.
           } # Placeholder fix as I reconstruct
           
           for (k in 1:K) {
             me_samples_k <- rowMeans(ey_plus[,,k] - ey_base[,,k], na.rm = TRUE)
             results_list[[length(results_list) + 1]] <- data.frame(
               Response         = resp,
               Predictor        = f_var,
               Category         = resp_levels[k],
               Effect           = mean(me_samples_k, na.rm = TRUE),
               Lower            = quantile(me_samples_k, (1 - prob) / 2, na.rm = TRUE),
               Upper            = quantile(me_samples_k, 1 - (1 - prob) / 2, na.rm = TRUE),
               Family           = dist,
               stringsAsFactors = FALSE
             )
           }
         } else {
           me_samples <- rowMeans(ey_plus - ey_base, na.rm = TRUE)
           results_list[[length(results_list) + 1]] <- data.frame(
             Response         = resp,
             Predictor        = f_var,
             Category         = NA,
             Effect           = mean(me_samples, na.rm = TRUE),
             Lower            = quantile(me_samples, (1 - prob) / 2, na.rm = TRUE),
             Upper            = quantile(me_samples, 1 - (1 - prob) / 2, na.rm = TRUE),
             Family           = dist,
             stringsAsFactors = FALSE
           )
         }
       }
    }
  }

  res_df <- do.call(rbind, results_list)
  if (!is.null(res_df)) {
     res_df$Response <- as.character(res_df$Response)
     res_df$Predictor <- as.character(res_df$Predictor)
     res_df$Family <- as.character(res_df$Family)
  }
  rownames(res_df) <- NULL
  return(res_df)
}
