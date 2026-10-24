#' Ensure Log-Likelihoods are Available
#'
#' Internal helper to check if a model object has pointwise log-likelihoods.
#' If not, it automatically refits the model (with reduced settings) to calculate them.
#'
#' @param model A fitted \code{because} model object.
#' @param quiet Logical; if \code{TRUE}, suppresses messages.
#'
#' @return A matrix of pointwise log-likelihoods (samples x observations).
#' @keywords internal
ensure_log_lik <- function(model, quiet = FALSE) {
    # 1. Check if log_lik exists in current samples
    if (!is.null(model$samples)) {
        col_names <- colnames(model$samples[[1]])
        log_lik_cols <- grep("^log_lik", col_names, value = TRUE)

        if (length(log_lik_cols) > 0) {
            # Found! Extract and return
            log_lik <- do.call(
                rbind,
                lapply(model$samples, function(chain) {
                    chain[, log_lik_cols, drop = FALSE]
                })
            )
            rownames(log_lik) <- NULL
            return(log_lik)
        }
    }

    # Pointwise log-likelihoods not found in current samples.
    if (!quiet) {
        message("Pointwise log-likelihoods not found in model samples.")
    }

    # Determine monitor names first (needed for both paths)
    monitor_names <- "log_lik" # Default fallback

    if (!is.null(model$parameter_map)) {
        responses <- unique(model$parameter_map$response)
        if (length(responses) > 0) {
            monitor_names <- paste0("log_lik_", responses)
        }
    }

    # 3. Attempt to use existing JAGS model object (Faster)
    # This works if the model was run sequentially and the session is active.
    # It fails if run in parallel (pointers often invalid) or loaded from disk.
    if (!is.null(model$model) && inherits(model$model, "jags")) {
        tryCatch(
            {
                if (!quiet) {
                    message(
                        "Attempting to use existing model state to calculate 'log_lik'..."
                    )
                }

                # Use coda.samples to extend the chain
                # Note: This advances the chain from where it left off.
                samples <- rjags::coda.samples(
                    model$model,
                    variable.names = monitor_names,
                    n.iter = 1000, # Sufficient for WAIC/LOO
                    thin = 1
                )

                if (!quiet) {
                    message("Calculated 'log_lik' from existing model.")
                }

                log_lik <- as.matrix(samples[[1]])
                rownames(log_lik) <- NULL
                return(log_lik)
            },
            error = function(e) {
                # Fallback to refit if any error occurs (e.g. dead pointer)
                if (!quiet) {
                    message(
                        "Existing model object invalid or unusable. Falling back to refit."
                    )
                }
            }
        )
    }

    # 4. Fallback: Full Refit (Recompile)
    if (!quiet) {
        message(
            "Refitting model (short chain) to calculate 'log_lik' for WAIC/LOO..."
        )
    }

    if (is.null(model$data) || is.null(model$model_code)) {
        stop(
            "Cannot refit model: 'data' or 'model_code' missing from model object."
        )
    }

    model_string <- paste(model$model_code, collapse = "\n")

    # Setup minimal run for efficiency
    n_chains <- 1
    n_adapt <- 500
    n_burnin <- 500
    n_iter <- 2000 # Results in 1500 samples
    n_thin <- 1

    # Compile
    jags_model <- rjags::jags.model(
        file = textConnection(model_string),
        data = model$data,
        n.chains = n_chains,
        n.adapt = n_adapt,
        quiet = TRUE
    )

    # Burn-in
    stats::update(jags_model, n.iter = n_burnin)

    # Sample
    samples <- rjags::coda.samples(
        jags_model,
        variable.names = monitor_names,
        n.iter = n_iter - n_burnin,
        thin = n_thin
    )

    if (!quiet) {
        message("Refitting complete.")
    }

    # Extract matrix
    log_lik <- as.matrix(samples[[1]])
    rownames(log_lik) <- NULL

    return(log_lik)
}
