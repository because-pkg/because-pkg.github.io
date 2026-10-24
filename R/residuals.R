#' @title Extract Fitted Values and Residuals for because Models
#' @name residuals.because
#' @description
#' Purely Bayesian methods for extracting expected fitted values and residuals
#' from fitted \code{because} model objects.
#'
#' Supports:
#' \itemize{
#'   \item \code{"response"}: Raw posterior residuals \eqn{e^{(s)} = y - \mu^{(s)}}
#'   \item \code{"pearson"}: Standardized residuals \eqn{(y - \mu^{(s)}) / \sigma^{(s)}}
#'   \item \code{"quantile"}: Randomized quantile residuals (Dunn & Smyth 1996; Hartig 2022)
#'         via the Probability Integral Transform (PIT) against \eqn{y^{\text{rep}}}.
#' }
#'
#' @param object A \code{because} model fit object.
#' @param resp Character string; name of the response variable. If \code{NULL}, takes the first response.
#' @param newdata Optional new dataset (currently evaluated within-sample).
#' @param ndraws Integer; number of posterior draws to use. Defaults to all draws (or 250 for quantile residuals).
#' @param re_formula Formula or \code{NA}; determines whether random/phylogenetic effects are included.
#'   If \code{NULL} (default), conditional on random effects.
#'   If \code{NA}, marginal prediction without random effects.
#' @param summary Logical; if \code{TRUE} (default), returns a summarized vector (posterior mean or median).
#'   If \code{FALSE}, returns the full \code{[ndraws x N_obs]} posterior matrix.
#' @param type Character; type of residual: \code{"quantile"} (default), \code{"response"}, or \code{"pearson"}.
#' @param ... Additional arguments.
#'
#' @references
#' Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
#'   \emph{Journal of Computational and Graphical Statistics}, 5(3), 236-244.
#'
#' Hartig, F. (2022). DHARMa: Residual Diagnostics for Hierarchical Regression Models.
#'
#' Gelman, A., & Hill, J. (2007). \emph{Data Analysis Using Regression and Multilevel/Hierarchical Models}.
#'   Cambridge University Press.
#'
#' @importFrom stats qnorm runif median
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' @export
fitted.because <- function(object, resp = NULL, newdata = NULL, ndraws = NULL, re_formula = NULL, summary = TRUE, ...) {
  # 1. Identification
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  # 2. Extract Posterior Samples
  samples_mat <- as.matrix(object$samples)
  n_total <- nrow(samples_mat)
  if (!is.null(ndraws) && ndraws < n_total) {
    set.seed(42)
    idx <- round(seq(1, n_total, length.out = ndraws))
    samples_mat <- samples_mat[idx, , drop = FALSE]
  }
  n_s <- nrow(samples_mat)
  
  # 3. Model Metadata & Resolution
  pm <- object$parameter_map
  pm_resp <- pm[pm$response == resp, ]
  if (nrow(pm_resp) == 0) {
    stop(paste("Response variable", resp, "not found in model parameters."))
  }
  
  data_list <- object$data
  obs_y <- if (!is.null(newdata)) {
    if (is.data.frame(newdata)) newdata[[resp]] else if (is.list(newdata) && !is.data.frame(newdata)) {
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
    object$original_data[[resp]] %||% data_list[[resp]]
  }
  
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
    search_data <- if (!is.null(newdata)) newdata else object$original_data
    if (idx_var %in% names(search_data)) {
      n_obs <- length(unique(search_data[[idx_var]]))
    } else if (paste0(idx_var, "_idx") %in% names(data_list)) {
      n_obs <- max(data_list[[paste0(idx_var, "_idx")]], na.rm = TRUE)
    }
  }
  
  # 4. Reconstruct Linear Predictor (eta)
  eta <- matrix(0, nrow = n_s, ncol = n_obs)
  
  alpha_param <- pm_resp$parameter[pm_resp$predictor == "(Intercept)"]
  if (length(alpha_param) > 0 && alpha_param %in% colnames(samples_mat)) {
    eta <- eta + replicate(n_obs, samples_mat[, alpha_param])
  }
  
  beta_rows <- pm_resp[pm_resp$type == "coefficient" & pm_resp$predictor != "(Intercept)", ]
  for (i in seq_len(nrow(beta_rows))) {
    p_name <- beta_rows$parameter[i]
    v_name <- beta_rows$predictor[i]
    if (p_name %in% colnames(samples_mat)) {
      b_samples <- samples_mat[, p_name]
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
      
      if (!is.null(x_vals)) {
        len_x <- if (is.matrix(x_vals)) ncol(x_vals) else length(x_vals)
        if (len_x > n_obs && target_level != "obs") {
          idx_var <- h_info$link_vars[[target_level]]
          links <- NULL
          for (df_name in names(search_data)) {
            if (is.data.frame(search_data[[df_name]]) && idx_var %in% names(search_data[[df_name]])) {
              links <- search_data[[df_name]][[idx_var]]
              if (length(links) == len_x) break else links <- NULL
            }
          }
          if (!is.null(links)) {
            first_idx <- !duplicated(links)
            x_vals <- x_vals[first_idx]
          } else {
            x_vals <- x_vals[1:n_obs]
          }
        }
        
        contribution <- if (is.matrix(x_vals)) {
          (as.numeric(b_samples) * x_vals)
        } else {
          (b_samples %*% t(as.matrix(x_vals)))
        }
        eta <- eta + contribution
      }
    }
  }
  
  # Random effects (if conditional)
  if (is.null(re_formula) || !is.na(re_formula)) {
    u_prefix <- paste0("u_", resp, "_")
    u_cols <- grep(paste0("^", u_prefix), colnames(samples_mat), value = TRUE)
    if (length(u_cols) > 0) {
      u_base_names <- unique(gsub("\\[\\d+\\]$", "", u_cols))
      for (ub in u_base_names) {
        lvl_name <- sub(u_prefix, "", ub)
        if (lvl_name == target_level || (lvl_name == "phylo" && target_level == "species")) {
          indices <- 1:n_obs
        } else {
          idx_name <- NULL
          if (!is.null(h_info) && !is.null(h_info$link_vars)) {
            if (lvl_name %in% names(h_info$link_vars)) {
              idx_name <- paste0(h_info$link_vars[[lvl_name]], "_idx")
            } else if (lvl_name == "phylo" && "species" %in% names(h_info$link_vars)) {
              idx_name <- paste0(h_info$link_vars[["species"]], "_idx")
            }
          }
          if (is.null(idx_name) || !idx_name %in% names(data_list)) {
            idx_name <- paste0(lvl_name, "_idx")
          }
          indices <- if (idx_name %in% names(data_list)) data_list[[idx_name]] else NULL
        }
        
        if (!is.null(indices)) {
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
  }
  
  # 5. Inverse Link Function
  fam <- if (!is.null(object$family) && resp %in% names(object$family)) object$family[[resp]] else "gaussian"
  mu_mat <- switch(fam,
    "gaussian"    = eta,
    "poisson"     = exp(eta),
    "binomial"    = 1 / (1 + exp(-eta)),
    "negbinomial" = exp(eta),
    "gamma"       = exp(eta),
    "lognormal"   = exp(eta),
    eta
  )
  
  if (summary) {
    return(colMeans(mu_mat, na.rm = TRUE))
  } else {
    return(mu_mat)
  }
}

#' @rdname residuals.because
#' @export
residuals.because <- function(object, resp = NULL, type = c("quantile", "response", "pearson"),
                              summary = TRUE, ndraws = 250, re_formula = NULL, ...) {
  type <- match.arg(type)
  
  # 1. Identification
  if (is.null(resp)) {
    resp <- as.character(all.vars(object$equations[[1]][[2]])[1])
  }
  
  # 2. Extract Observed Data
  y <- object$original_data[[resp]] %||% object$data[[resp]]
  if (is.null(y)) {
    stop(paste("Could not find observed data for response:", resp))
  }
  
  # Handle hierarchical resolution alignment
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
  
  # 3. Compute Residuals Based on Type
  if (type == "quantile") {
    # Dunn & Smyth (1996) Randomized Quantile Residuals via Posterior Predictive Draws
    yrep <- posterior_predict(object, resp = resp, ndraws = ndraws, re_formula = re_formula)
    
    # Align dimensions if y is observation-level but yrep is higher level
    if (ncol(yrep) != length(y)) {
      if (length(y) > ncol(yrep) && target_level != "obs" && !is.null(h_info$link_vars[[target_level]])) {
        idx_var <- h_info$link_vars[[target_level]]
        search_data <- object$original_data
        links <- if (idx_var %in% names(search_data)) search_data[[idx_var]] else object$data[[paste0(idx_var, "_idx")]]
        if (!is.null(links) && length(links) == length(y)) {
          y <- y[!duplicated(links)]
        } else {
          y <- y[seq_len(ncol(yrep))]
        }
      } else {
        yrep <- yrep[, seq_along(y), drop = FALSE]
      }
    }
    
    S <- nrow(yrep)
    N <- length(y)
    u <- numeric(N)
    
    for (i in seq_len(N)) {
      p_lower <- mean(yrep[, i] < y[i], na.rm = TRUE)
      p_equal <- mean(yrep[, i] == y[i], na.rm = TRUE)
      
      # Probability Integral Transform with randomized jitter for discrete ties
      u_val <- p_lower + stats::runif(1) * p_equal
      # Boundary protection to avoid infinite qnorm values
      u_val <- min(max(u_val, 1 / (2 * S)), 1 - 1 / (2 * S))
      u[i] <- u_val
    }
    
    # Map to standard normal scale
    res_val <- stats::qnorm(u)
    return(res_val)
    
  } else {
    # Continuous / Response or Pearson Residuals
    mu_mat <- fitted.because(object, resp = resp, ndraws = ndraws, re_formula = re_formula, summary = FALSE)
    
    if (ncol(mu_mat) != length(y)) {
      if (length(y) > ncol(mu_mat)) {
        y <- y[seq_len(ncol(mu_mat))]
      } else {
        mu_mat <- mu_mat[, seq_along(y), drop = FALSE]
      }
    }
    
    # Raw response residuals matrix: y - mu
    y_mat <- matrix(rep(y, each = nrow(mu_mat)), nrow = nrow(mu_mat))
    res_mat <- y_mat - mu_mat
    
    if (type == "pearson") {
      # Standardize by posterior sigma / dispersion
      samples_mat <- as.matrix(object$samples)
      if (nrow(samples_mat) > nrow(mu_mat)) {
        set.seed(42)
        idx <- round(seq(1, nrow(samples_mat), length.out = nrow(mu_mat)))
        samples_mat <- samples_mat[idx, , drop = FALSE]
      }
      
      sigma_name <- paste0("sigma_", resp, "_res")
      tau_name <- paste0("tau_res_", resp)
      sigma_s <- if (sigma_name %in% colnames(samples_mat)) {
        samples_mat[, sigma_name]
      } else if (tau_name %in% colnames(samples_mat)) {
        1 / sqrt(samples_mat[, tau_name])
      } else {
        rep(1, nrow(mu_mat))
      }
      
      res_mat <- res_mat / matrix(rep(sigma_s, ncol(res_mat)), ncol = ncol(res_mat))
    }
    
    if (summary) {
      return(colMeans(res_mat, na.rm = TRUE))
    } else {
      return(res_mat)
    }
  }
}
