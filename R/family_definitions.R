#' @title Family Definition Generics and Default Methods
#' @description S3 generics and methods for distribution families in because.
#' This enables custom distributions to be added by defining S3 methods.
#' @name family_definitions
#' @importFrom stats rpois rbinom runif rnbinom rnorm dpois
# --- Global NIMBLE Distributions for because package ---

# These are defined lazily to avoid a hard dependency on nimble during package installation.
.nimble_fn_cache <- new.env(parent = emptyenv())

#' @keywords internal
get_nimble_fn <- function(name) {
    if (!requireNamespace("nimble", quietly = TRUE)) {
        stop("Package 'nimble' is required for this operation. Please install it with install.packages('nimble').")
    }
    
    if (exists(name, envir = .nimble_fn_cache)) {
        return(get(name, envir = .nimble_fn_cache))
    }
    
    fn <- switch(name,
        "dnb_because" = nimble::nimbleFunction(
            run = function(x = double(0), v_mu = double(0), v_r = double(0), log = integer(0)) {
                returnType(double(0))
                s_mu <- max(1e-10, min(1000000, v_mu))
                s_r <- max(1e-10, min(1000000, v_r))
                log_r_plus_mu <- log(s_r + s_mu)
                log_prob <- lgamma(x + s_r) - lgamma(x + 1) - lgamma(s_r) + 
                            s_r * log(s_r) + x * log(s_mu) - (x + s_r) * log_r_plus_mu
                if (log) return(log_prob) else return(exp(log_prob))
            }
        ),
        "rnb_because" = nimble::nimbleFunction(
            run = function(n = integer(0), v_mu = double(0), v_r = double(0)) {
                returnType(double(0))
                prob <- v_r / (v_r + v_mu)
                return(rnbinom(1, size = v_r, prob = prob))
            }
        ),
        "dzip_because" = nimble::nimbleFunction(
            run = function(x = double(0), v_mu = double(0), v_psi = double(0), log = integer(0)) {
                returnType(double(0))
                if (x == 0) {
                    val <- v_psi + (1 - v_psi) * exp(-v_mu)
                } else {
                    val <- (1 - v_psi) * dpois(x, v_mu, 0)
                }
                if (log) return(log(max(1.0e-30, val))) else return(val)
            }
        ),
        "rzip_because" = nimble::nimbleFunction(
            run = function(n = integer(0), v_mu = double(0), v_psi = double(0)) {
                returnType(double(0))
                if (runif(1) < v_psi) {
                    return(0)
                } else {
                    return(rpois(1, v_mu))
                }
            }
        ),
        "dzinb_because" = nimble::nimbleFunction(
            run = function(x = double(0), v_mu = double(0), v_r = double(0), v_psi = double(0), log = integer(0)) {
                returnType(double(0))
                s_mu <- max(1e-10, min(1000000, v_mu))
                s_r <- max(1e-10, min(1000000, v_r))
                if (x == 0) {
                    log_p_nb_zero <- s_r * (log(s_r) - log(s_r + s_mu))
                    val <- v_psi + (1 - v_psi) * exp(log_p_nb_zero)
                    log_prob <- log(max(1e-30, val))
                } else {
                    log_r_plus_mu <- log(s_r + s_mu)
                    log_nb <- lgamma(x + s_r) - lgamma(x + 1) - lgamma(s_r) + 
                                s_r * log(s_r) + x * log(s_mu) - (x + s_r) * log_r_plus_mu
                    log_prob <- log(max(1e-30, 1 - v_psi)) + log_nb
                }
                if (log) return(log_prob) else return(exp(log_prob))
            }
        ),
        "rzinb_because" = nimble::nimbleFunction(
            run = function(n = integer(0), v_mu = double(0), v_r = double(0), v_psi = double(0)) {
                returnType(double(0))
                if (runif(1) < v_psi) {
                    return(0)
                } else {
                    prob <- v_r / (v_r + v_mu)
                    return(rnbinom(1, size = v_r, prob = prob))
                }
            }
        )
    )
    
    assign(name, fn, envir = .nimble_fn_cache)
    return(fn)
}

# Dummy exports to satisfy NAMESPACE during installation/documentation
# In practice these should not be used directly by users without nimble.

#' @export
dnb_because <- function(...) { stop("Use nimble_family_optimization to access this NIMBLE function.") }
#' @export
rnb_because <- function(...) { stop("Use nimble_family_optimization to access this NIMBLE function.") }
#' @export
dzip_because <- function(...) { stop("Use nimble_family_optimization to access this NIMBLE function.") }
#' @export
rzip_because <- function(...) { stop("Use nimble_family_optimization to access this NIMBLE function.") }
#' @export
dzinb_because <- function(...) { stop("Use nimble_family_optimization to access this NIMBLE function.") }
#' @export
rzinb_because <- function(...) { stop("Use nimble_family_optimization to access this NIMBLE function.") }


# --- Standard family definitions ---
NULL

#' Get default precision prior for a distribution family
#'
#' @param family A family object
#' @param param_name Name of the parameter (e.g. "tau_e_Abundance")
#' @param ... Additional arguments
#' @return Character string; JAGS prior statement (e.g. "tau_e_Abundance ~ dgamma(1, 1)")
#' @keywords internal
#' @export
jags_family_precision_prior <- function(family, param_name, ...) {
    UseMethod("jags_family_precision_prior")
}

#' @keywords internal
#' @export
jags_family_precision_prior.default <- function(family, param_name, ...) {
    # Default to weakly informative but regularizing prior for non-Gaussian/unknown
    return(paste0(param_name, " ~ dgamma(10, 10)"))
}

#' @keywords internal
#' @export
jags_family_precision_prior.because_family_gaussian <- function(
    family,
    param_name,
    ...
) {
    # Half-uniform prior on sigma: ONLY for residual error parameters (tau_res_*).
    # This is scale-invariant and robust to raw or standardised data.
    # Gelman (2006) advises against priors directly on precision for residual variance.
    #
    # For all other Gaussian precision parameters (tau_u_* random effects,
    # tau_obs_* measurement error, tau_struct_* structural terms, etc.),
    # fall back to dgamma(1, 1) — appropriate weakly-informative regularisation
    # for variance components where the scale is already model-defined.
    if (grepl("^tau_res_", param_name)) {
        sigma_name <- sub("^tau_res_", "sigma_", param_name)
        sigma_name <- paste0(sigma_name, "_res")
        return(c(
            paste0(sigma_name, " ~ dunif(0, 100)"),
            paste0(param_name, " <- 1 / (", sigma_name, " * ", sigma_name, ")")
        ))
    }
    # Non-residual precision: weakly-informative Gamma prior
    return(paste0(param_name, " ~ dgamma(1, 1)"))
}


#' @keywords internal
#' @export
jags_family_precision_prior.because_family_poisson <- function(
    family,
    param_name,
    ...
) {
    # Regularizing prior for log-link stability
    return(paste0(param_name, " ~ dgamma(10, 10)"))
}

#' @keywords internal
#' @export
jags_family_precision_prior.because_family_binomial <- function(
    family,
    param_name,
    ...
) {
    # Regularizing prior for logit-link stability
    return(paste0(param_name, " ~ dgamma(10, 10)"))
}

#' Generate JAGS likelihood code for a distribution family
#'
#' This generic allows S3 dispatch to generate the appropriate JAGS
#' likelihood code for any distribution family.
#'
#' @param family A family object (created by \code{\link{get_family_object}} or custom constructor)
#' @param response Character string; name of the response variable
#' @param predictors Character vector; names of predictor variables (may be NULL)
#' @param suffix Character string; suffix for variable names (e.g., "1" for multiple responses)
#' @param has_structure Logical; whether the model includes a structure (e.g., phylogenetic)
#' @param link Character string; link function ("identity", "log", "logit")
#' @param ... Additional arguments passed to methods
#'
#' @return A list with:
#' \describe{
#'   \item{likelihood_code}{Character vector of JAGS likelihood statements}
#'   \item{prior_code}{Character vector of JAGS prior statements (or NULL)}
#'   \item{data_requirements}{Character vector of required data elements (or NULL)}
#' }
#'
#' @keywords internal
#' @export
jags_family_likelihood <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "identity",
    ...
) {
    UseMethod("jags_family_likelihood")
}

#' @keywords internal
#' @export
jags_family_likelihood.default <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "identity",
    ...
) {
    # Get family name
    fam_name <- if (is.character(family)) family else family$family
    stop(paste0(
        "No jags_family_likelihood method defined for family '",
        fam_name,
        "'.\n",
        "Use because_family() to create custom family methods, or check if ",
        "an extension module package is required."
    ))
}

#' @keywords internal
#' @export
jags_family_likelihood.because_family_gaussian <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "identity",
    ...
) {
    # Standard normal likelihood
    mu_var <- paste0("mu_", response, suffix)
    tau_var <- paste0("tau_res_", response, suffix)

    likelihood_code <- paste0(
        "    ",
        response,
        "[i] ~ dnorm(",
        mu_var,
        "[i], ",
        tau_var,
        ")"
    )

    prior_code <- jags_family_precision_prior(family, tau_var)

    list(
        likelihood_code = likelihood_code,
        prior_code = prior_code,
        data_requirements = NULL
    )
}

#' @keywords internal
#' @export
jags_family_likelihood.because_family_binomial <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "logit",
    ...
) {
    # Bernoulli likelihood with logit link
    p_var <- paste0("p_", response, suffix)
    mu_var <- paste0("mu_", response, suffix)

    likelihood_code <- c(
        paste0("    ", p_var, "[i] <- ilogit(", mu_var, "[i])"),
        paste0("    ", response, "[i] ~ dbern(", p_var, "[i])")
    )

    list(
        likelihood_code = likelihood_code,
        prior_code = NULL,
        data_requirements = NULL
    )
}

#' @keywords internal
#' @export
jags_family_likelihood.because_family_poisson <- function(
    family,
    response,
    predictors = NULL,
    suffix = "",
    has_structure = FALSE,
    link = "log",
    ...
) {
    # Poisson likelihood with log link
    lambda_var <- paste0("lambda_", response, suffix)
    mu_var <- paste0("mu_", response, suffix)

    likelihood_code <- c(
        paste0("    ", lambda_var, "[i] <- exp(", mu_var, "[i])"),
        paste0("    ", response, "[i] ~ dpois(", lambda_var, "[i])")
    )

    list(
        likelihood_code = likelihood_code,
        prior_code = NULL,
        data_requirements = NULL
    )
}

#' @keywords internal
#' @export
needs_zero_inflation_hook.because_family_zip <- function(
    family,
    variable_name,
    ...
) {
    return(TRUE)
}

#' @keywords internal
#' @export
needs_zero_inflation_hook.because_family_zinb <- function(
    family,
    variable_name,
    ...
) {
    return(TRUE)
}

#' @keywords internal
#' @export
nimble_family_optimization.because_family_zip <- function(
    family,
    model_string,
    variable,
    ...
) {
    list(
        model_string = model_string,
        nimble_functions = list(
            dzip_because = get_nimble_fn("dzip_because"),
            rzip_because = get_nimble_fn("rzip_because")
        )
    )
}

#' @keywords internal
#' @export
nimble_family_optimization.because_family_zinb <- function(
    family,
    model_string,
    variable,
    ...
) {
    list(
        model_string = model_string,
        nimble_functions = list(
            dzinb_because = get_nimble_fn("dzinb_because"),
            rzinb_because = get_nimble_fn("rzinb_because")
        )
    )
}

#' @keywords internal
#' @export
nimble_family_optimization.because_family_negbinomial <- function(
    family,
    model_string,
    variable,
    ...
) {
    list(
        model_string = model_string,
        nimble_functions = list(
            dnb_because = get_nimble_fn("dnb_because"),
            rnb_because = get_nimble_fn("rnb_because")
        )
    )
}

#' Create a family object for a given distribution
#'
#' Returns a family object that can be dispatched on via S3.
#'
#' @param name Character string; family name (e.g., "gaussian", "binomial")
#' @return A family object with appropriate S3 class
#' @examples
#' get_family("gaussian")
#' get_family("binomial")
#' @export
get_family <- function(name) {
    # Normalize name
    name <- tolower(name)

    # Map aliases
    name <- switch(name, "normal" = "gaussian", "bernoulli" = "binomial", name)

    # Create family object with S3 class
    structure(
        list(family = name),
        class = c(paste0("because_family_", name), "because_family")
    )
}
