#' Calculate WAIC with Standard Errors for a Because Model
#'
#' Calculates the Widely Applicable Information Criterion (WAIC) with standard errors
#' for a fitted Because model using pointwise log-likelihoods.
#'
#' @param model A fitted model object of class \code{"because"} returned by \code{\link{because}}
#'   with \code{WAIC = TRUE}.
#'
#' @return A data frame with columns \code{Estimate} and \code{SE} containing:
#' \item{elpd_waic}{Expected log pointwise predictive density (higher is better)}
#' \item{p_waic}{Effective number of parameters}
#' \item{waic}{The WAIC value (lower is better for model comparison)}
#'
#' The returned object also has a \code{pointwise} attribute containing individual
#' observation contributions for model comparison.
#'
#' @details
#' This function is automatically called by \code{\link{because}} when \code{WAIC = TRUE}.
#' It can also be called manually. If the model was not originally fitted with
#' \code{WAIC = TRUE} (so pointwise log-likelihoods are missing), this function
#' will automatically refit the model (using a short MCMC run) to compute them.
#'
#' @section WAIC Definition:
#' The Widely Applicable Information Criterion (WAIC) is calculated as:
#' \deqn{WAIC = -2 \times (lppd - p_{waic})}
#' where:
#' \itemize{
#'   \item \eqn{lppd = \sum_{i=1}^N \log(\frac{1}{S} \sum_{s=1}^S \exp(log\_lik_{is}))} is the log pointwise predictive density
#'   \item \eqn{p_{waic} = \sum_{i=1}^N \text{var}(log\_lik_{is})} is the effective number of parameters
#' }
#'
#' **WAIC Algorithm**:
#' \enumerate{
#'   \item \strong{lpd} (log pointwise predictive density): For each observation \eqn{i},
#'     compute \eqn{\log(\text{mean}(\exp(\text{log\_lik}_i)))} across MCMC samples
#'   \item \strong{p_waic}: For each observation \eqn{i}, compute
#'     \eqn{\text{var}(\text{log\_lik}_i)} across MCMC samples
#'   \item \strong{elpd_waic}: \eqn{\text{lpd}_i - \text{p\_waic}_i} for each observation
#'   \item \strong{waic}: \eqn{-2 \times \sum \text{elpd\_waic}_i}
#' }
#'
#' @section Standard Errors:
#' Standard errors for WAIC are calculated using the pointwise contributions:
#' \deqn{SE(WAIC) = \sqrt{N \times \text{var}(waic_i)}}
#' where \eqn{waic_i = -2 \times (lppd_i - p_{waic,i})}.
#'
#' @examples
#' \dontrun{
#'   # Fit model with WAIC monitoring
#'   fit <- because(data, tree, equations, WAIC = TRUE)
#'
#'   # View WAIC with standard errors
#'   fit$WAIC
#'   #             Estimate   SE
#'   # elpd_waic   -617.3   12.4
#'   # p_waic        12.3    3.1
#'   # waic        1234.5   24.8
#'
#'   # Compare two models
#'   fit1$WAIC
#'   fit2$WAIC
#'   # Model with lower WAIC is preferred
#'   # Difference is significant if |WAIC1 - WAIC2| > 2 * sqrt(SE1^2 + SE2^2)
#' }
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC.
#' \emph{Statistics and Computing}, 27(5), 1413-1432.
#'
#' @export
because_waic <- function(model) {
    if (!inherits(model, "because")) {
        stop("Input must be a 'because' model object.")
    }

    # Extract log_lik using helper (auto-refits if needed)
    log_lik <- ensure_log_lik(model)

    # Number of observations and samples
    n_obs <- ncol(log_lik)
    n_samples <- nrow(log_lik)

    if (n_samples < 100) {
        warning(
            "Only ",
            n_samples,
            " MCMC samples available. ",
            "Standard errors may be unreliable. Consider using more iterations."
        )
    }

    # Compute WAIC following Vehtari et al. (2017)

    # 1. Compute lpd (log pointwise predictive density)
    #    For each observation: log(mean(exp(log_lik)))
    lpd_i <- apply(log_lik, 2, function(ll_i) {
        # Use log-sum-exp trick for numerical stability
        max_ll <- max(ll_i)
        log(mean(exp(ll_i - max_ll))) + max_ll
    })

    # 2. Compute p_waic (effective number of parameters)
    #    For each observation: var(log_lik)
    p_waic_i <- apply(log_lik, 2, var)

    # 3. Compute elpd_waic (expected log pointwise predictive density)
    elpd_waic_i <- lpd_i - p_waic_i

    # 4. Compute pointwise WAIC
    waic_i <- -2 * elpd_waic_i

    # 5. Sum across observations for totals
    elpd_waic <- sum(elpd_waic_i)
    p_waic <- sum(p_waic_i)
    waic <- sum(waic_i)

    # 6. Compute standard errors
    #    SE = sqrt(N * var(pointwise contributions))
    se_elpd_waic <- sqrt(n_obs * var(elpd_waic_i))
    se_p_waic <- sqrt(n_obs * var(p_waic_i))
    se_waic <- sqrt(n_obs * var(waic_i))

    # 7. Create result data frame
    result <- data.frame(
        Estimate = c(elpd_waic, p_waic, waic),
        SE = c(se_elpd_waic, se_p_waic, se_waic),
        row.names = c("elpd_waic", "p_waic", "waic")
    )

    # 8. Store pointwise values as attribute for model comparison
    attr(result, "pointwise") <- data.frame(
        elpd_waic = elpd_waic_i,
        p_waic = p_waic_i,
        waic = waic_i
    )

    # 9. Store dimensions
    attr(result, "dims") <- c(n_obs = n_obs, n_samples = n_samples)

    class(result) <- c("because_waic", "data.frame")

    return(result)
}

#' @keywords internal
#' @export
print.because_waic <- function(x, digits = 1, ...) {
    cat("WAIC with Standard Errors\n")
    cat("-------------------------\n")

    dims <- attr(x, "dims")
    if (!is.null(dims)) {
        cat(sprintf(
            "N = %d observations, %d MCMC samples\n\n",
            dims["n_obs"],
            dims["n_samples"]
        ))
    }

    # Print as data frame with custom formatting
    x_print <- x
    x_print$Estimate <- round(x_print$Estimate, digits)
    x_print$SE <- round(x_print$SE, digits)

    print.data.frame(x_print, ...)

    cat("\nInterpretation:\n")
    cat("  - Lower WAIC indicates better model fit\n")
    cat(
        "  - Compare models: difference significant if |WAIC1 - WAIC2| > 2*SE_diff\n"
    )
    cat("  - p_waic estimates effective number of parameters\n")

    invisible(x)
}
