#' Causal Interventions for Structural Equation Models
#'
#' Evaluates counterfactuals using Pearl's do-operator. Supports atomic
#' interventions (setting a variable to a fixed value), shift interventions
#' (adding a constant to the natural value), and stochastic interventions
#' (drawing from a distribution).
#'
#' @param object A fitted model object.
#' @param ... Interventions. Can be numeric constants, formulas, or expressions.
#' @return Intervened data or counterfactual samples.
#' @export
do <- function(object, ...) {
    UseMethod("do")
}

#' @rdname do
#' @export
because_do <- do

#' Causal Interventions for Because Models
#'
#' This method implements Pearl's do-operator for Bayesian Structural Equation 
#' Models fitted with \code{because}. It topologically propagates interventions 
#' downstream through the causal graph to simulate counterfactual outcomes.
#'
#' @param object A fitted \code{because} model object.
#' @param ... Interventions specified as name-value pairs. Supported types include:
#'   \itemize{
#'     \item \strong{Atomic (Hard)}: \code{temp = 25} (fixes value exactly)
#'     \item \strong{Shift (Additive)}: \code{temp = ~ . + 1.5} (adds 1.5 to natural values)
#'     \item \strong{Percentage Shift}: \code{temp = "+5\%"} or \code{temp = "-10\%"} (multiplies by 1.05 or 0.90)
#'     \item \strong{Stochastic}: \code{temp = ~ rnorm(n, mean = . + 1.5, sd = 0.05)}
#'   }
#'   Note: In formulas, \code{.} evaluates to the current/natural values, and \code{n} 
#'   evaluates to the number of elements required for random generator functions.
#' @param ndraws Integer; number of posterior draws to simulate. Defaults to all draws.
#' @param re_formula Formula; determines which random effects to condition on.
#'   Defaults to \code{NULL} (conditional on estimated groups, e.g., historical sites).
#'   Set to \code{NA} for marginal predictions over new, unmeasured groups.
#' @param raw_scale Logical; if \code{TRUE}, the engine will automatically detect variables 
#'   that were z-scored using \code{scale()} prior to model fitting. It will unscale the 
#'   data, apply your intervention on the raw metric, and then re-scale it back to the 
#'   z-metric before propagating the counterfactual. Defaults to \code{FALSE}.
#'   \strong{Important Note}: Base R's \code{data.frame(X = scale(X))} will silently strip 
#'   the scaling attributes required for this feature to work. To preserve them, you must 
#'   assign the scaled column to an existing data frame, e.g., \code{df$X <- scale(X)}.
#' @return An object of class \code{because_counterfactual} (a list) containing 
#'   \code{[ndraws x N_obs]} matrices for each endogenous variable under the intervention.
#'
#' @section Warning on Scaled Data:
#'   If you scaled your data (e.g., z-scoring) prior to fitting the model to obtain 
#'   comparable path coefficients, **interventions must be specified on the scaled metric** 
#'   unless you set \code{raw_scale = TRUE}. For example, if Temperature was z-scored, 
#'   \code{do(fit, Temp = ~ . + 1.5, raw_scale = FALSE)} increases Temperature by 
#'   **1.5 standard deviations**, not 1.5 degrees. 
#'   
#'   Furthermore, applying **percentage shifts** (e.g., \code{Temp = "+5\%"}) to z-scored 
#'   data will multiply the z-scores by 1.05. Because the mean of z-scored data is 0, 
#'   this stretches the variance but does *not* shift the mean (since 0 * 1.05 = 0). 
#'   Always carefully consider the scale of your variables when designing interventions, 
#'   and utilize \code{raw_scale = TRUE} if you wish to intervene on the original metric.
#'
#' @examples
#' \dontrun{
#' # Fit a simple mediation model
#' df <- data.frame(X = rnorm(100), M = rnorm(100), Y = rnorm(100))
#' fit <- because(list(M ~ X, Y ~ M + X), data = df)
#'
#' # 1. Hard intervention: set X to exactly 15 everywhere
#' res_hard <- do(fit, X = 15)
#' mean(res_hard$Y)
#'
#' # 2. Shift intervention: increase X by 5 from its baseline
#' res_shift <- do(fit, X = ~ . + 5)
#'
#' # 3. Stochastic shift: increase X by 5, with variance
#' res_stoch <- do(fit, X = ~ rnorm(n, mean = . + 5, sd = 0.5))
#' }
#' @export
do.because <- function(object, ..., ndraws = NULL, re_formula = NULL, raw_scale = FALSE) {
    interventions <- list(...)
    if (length(interventions) == 0) {
        warning("No interventions specified. Returning standard predictions.")
    }
    
    # Extract total draws
    n_total <- nrow(as.matrix(object$samples))
    if (is.null(ndraws)) ndraws <- n_total
    
    # 1. Get equations and variables
    eqs <- object$equations
    # Build simple adjacency list
    adj <- list()
    all_vars <- c()
    for (eq in eqs) {
        lhs <- as.character(eq[[2]])
        rhs <- all.vars(eq[[3]])
        adj[[lhs]] <- rhs
        all_vars <- unique(c(all_vars, lhs, rhs))
    }
    
    # 2. Topological Sort
    topo_order <- c()
    visited <- setNames(rep(FALSE, length(all_vars)), all_vars)
    temp_mark <- setNames(rep(FALSE, length(all_vars)), all_vars)
    
    visit <- function(n) {
        if (temp_mark[[n]]) stop("Cyclic graph detected.")
        if (!visited[[n]]) {
            temp_mark[[n]] <<- TRUE
            # Visit parents
            parents <- adj[[n]]
            for (p in parents) {
                if (!is.null(p) && p %in% names(visited)) visit(p)
            }
            temp_mark[[n]] <<- FALSE
            visited[[n]] <<- TRUE
            topo_order <<- c(topo_order, n)
        }
    }
    for (v in all_vars) visit(v)
    
    # 3. Initialize newdata with original observed data
    # We expand everything to [ndraws x n_obs] matrices as we go, or keep them as vectors
    # until they are modified.
    newdata <- list()
    
    # Need to handle original_data structure which might be hierarchical
    # Flatten it into a single list for easy lookup, but preserve row counts
    flat_data <- list()
    if (is.data.frame(object$original_data)) {
        flat_data <- as.list(object$original_data)
    } else if (is.list(object$original_data)) {
        for (df_name in names(object$original_data)) {
            if (is.data.frame(object$original_data[[df_name]])) {
                for (cn in names(object$original_data[[df_name]])) {
                    if (is.null(flat_data[[cn]])) {
                        flat_data[[cn]] <- object$original_data[[df_name]][[cn]]
                    }
                }
            } else {
                flat_data[[df_name]] <- object$original_data[[df_name]]
            }
        }
    }
    
    newdata <- flat_data
    
    # 4. Simulation Loop
    results <- list()
    
    for (var in topo_order) {
        n_obs <- length(newdata[[var]])
        if (is.matrix(newdata[[var]])) n_obs <- ncol(newdata[[var]])
        if (is.null(n_obs)) next # Skip latent or unmeasured root nodes
        
        # Is there an intervention on this variable?
        if (var %in% names(interventions)) {
            iv <- interventions[[var]]
            current_val <- newdata[[var]]
            
            # --- Apply raw_scale Unscaling if requested ---
            unscaled <- FALSE
            if (raw_scale && !is.null(object$scale_info[[var]])) {
                s_info <- object$scale_info[[var]]
                current_val <- (current_val * s_info$scale) + s_info$center
                unscaled <- TRUE
            }
            
            if (is.numeric(iv) && length(iv) == 1) {
                # Hard / Atomic intervention
                sim_val <- matrix(iv, nrow = ndraws, ncol = n_obs)
            } else if (is.character(iv) && length(iv) == 1 && grepl("%$", trimws(iv))) {
                # Percentage shift intervention (e.g. "+5%", "-10%")
                val_str <- trimws(iv)
                val_num <- as.numeric(sub("%$", "", val_str))
                if (is.na(val_num)) stop(paste("Invalid percentage format for", var, ":", iv))
                
                multiplier <- 1 + (val_num / 100)
                
                base_mat <- current_val
                if (!is.matrix(base_mat)) {
                    base_mat <- matrix(rep(base_mat, each = ndraws), nrow = ndraws, ncol = n_obs)
                }
                sim_val <- base_mat * multiplier
            } else if (inherits(iv, "formula")) {
                # Shift or Stochastic intervention via formula (e.g. ~ . + 1.5)
                # Convert . to current value
                expr <- iv[[2]]
                
                # If current_val is a vector, we need to broadcast it to a matrix
                # so the formula applies row-wise across draws
                base_mat <- current_val
                if (!is.matrix(base_mat)) {
                    base_mat <- matrix(rep(base_mat, each = ndraws), nrow = ndraws, ncol = n_obs)
                }
                
                # Evaluate expression
                # We provide '.' as the base matrix, and 'n' as length for rnorm
                eval_env <- list(. = base_mat, n = length(base_mat))
                
                # Expose all other upstream variables to the formula just in case
                for (vname in names(newdata)) {
                    eval_env[[vname]] <- newdata[[vname]]
                }
                
                sim_val <- eval(expr, envir = eval_env)
                
                # If the function returned a vector of length n_obs (like just a constant shift),
                # or length ndraws * n_obs, ensure it's a matrix
                if (!is.matrix(sim_val)) {
                    if (length(sim_val) == n_obs) {
                         sim_val <- matrix(rep(sim_val, each = ndraws), nrow = ndraws, ncol = n_obs)
                    } else if (length(sim_val) == ndraws * n_obs) {
                         sim_val <- matrix(sim_val, nrow = ndraws, ncol = n_obs)
                    } else {
                         stop(paste("Intervention formula for", var, "returned unexpected dimensions."))
                    }
                }
            } else {
                stop(paste("Unsupported intervention type for", var))
            }
            
            # --- Apply raw_scale Rescaling ---
            if (unscaled) {
                sim_val <- (sim_val - s_info$center) / s_info$scale
            }
            
            newdata[[var]] <- sim_val
            results[[var]] <- sim_val
            
        } else if (var %in% names(adj)) {
            # Not intervened, but it is an endogenous variable (has an equation)
            # Predict it using current upstream newdata
            pred <- posterior_predict(object, resp = var, newdata = newdata, ndraws = ndraws, re_formula = re_formula)
            newdata[[var]] <- pred
            results[[var]] <- pred
        }
        # If it's neither intervened nor endogenous, it's a root node, keep original values.
    }
    
    # --- Final Unscaling ---
    if (raw_scale) {
        for (var in names(results)) {
            if (!is.null(object$scale_info[[var]])) {
                s_info <- object$scale_info[[var]]
                results[[var]] <- (results[[var]] * s_info$scale) + s_info$center
            }
        }
    }
    
    class(results) <- c("because_counterfactual", "list")
    return(results)
}

#' Summarize a Because Counterfactual Object
#'
#' @param object An object of class \code{because_counterfactual} returned by \code{do()}.
#' @param ... Additional arguments (ignored).
#' @return A data.frame summarizing the expected value and credible intervals for each variable.
#' @export
summary.because_counterfactual <- function(object, ...) {
    res <- data.frame()
    
    for (var in names(object)) {
        mat <- object[[var]]
        if (is.matrix(mat)) {
            # Compute global average per posterior draw first
            # (i.e. average across sites/observations for each draw)
            global_draws <- rowMeans(mat)
            
            row_df <- data.frame(
                Variable = var,
                Mean = mean(global_draws),
                SD = sd(global_draws),
                `Q2.5` = stats::quantile(global_draws, 0.025, names = FALSE),
                `Q50` = stats::quantile(global_draws, 0.50, names = FALSE),
                `Q97.5` = stats::quantile(global_draws, 0.975, names = FALSE),
                stringsAsFactors = FALSE
            )
            # Fix column names that might get mangled by data.frame
            names(row_df) <- c("Variable", "Mean", "SD", "2.5%", "50%", "97.5%")
            res <- rbind(res, row_df)
        }
    }
    rownames(res) <- NULL
    class(res) <- c("summary_because_counterfactual", "data.frame")
    return(res)
}

#' @export
print.summary_because_counterfactual <- function(x, ...) {
    cat("Counterfactual Simulation Summary\n")
    cat("---------------------------------\n")
    cat("Estimates represent the global expectation (averaged across all observations)\n")
    cat("under the intervened causal structure.\n\n")
    
    class(x) <- "data.frame"
    print(x, row.names = FALSE, digits = 4)
    invisible(x)
}
