#' @title Family Definition Generics and Default Methods
#' @description S3 generics and methods for distribution families in because.
#' This enables custom distributions to be added by defining S3 methods.
#' @importFrom stats rpois rbinom runif
#' @importFrom nimble nimbleFunction returnType
NULL

# --- Global NIMBLE Distributions for because package ---

#' @export
dnb_because <- nimble::nimbleFunction(
    run = function(
        x = double(0),
        mu = double(0),
        r = double(0),
        log = integer(0, default = 0)
    ) {
        returnType(double(0))
        # Use R-style dnbinom which is known to the NIMBLE compiler
        prob <- dnbinom(x, size = r, mu = mu, log = 0)
        if (log) return(log(max(1.0e-30, prob))) else return(prob)
    }
)

#' @export
rnb_because <- nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), r = double(0)) {
        returnType(double(0))
        # Use R-style rnbinom which is known to the NIMBLE compiler
        return(rnbinom(1, size = r, mu = mu))
    }
)

#' @export
dzip_because <- nimble::nimbleFunction(
    run = function(
        x = double(0),
        mu = double(0),
        psi = double(0),
        log = integer(0, default = 0)
    ) {
        returnType(double(0))
        if (x == 0) {
            prob <- psi + (1 - psi) * exp(-mu)
        } else {
            prob <- (1 - psi) * dpois(x, mu, 0)
        }
        if (log) return(log(max(1.0e-30, prob))) else return(prob)
    }
)

#' @export
rzip_because <- nimble::nimbleFunction(
    run = function(n = integer(0), mu = double(0), psi = double(0)) {
        returnType(double(0))
        if (runif(1) < psi) {
            return(0)
        } else {
            return(rpois(1, mu))
        }
    }
)

#' @export
dzinb_because <- nimble::nimbleFunction(
    run = function(
        x = double(0),
        mu = double(0),
        r = double(0),
        psi = double(0),
        log = integer(0, default = 0)
    ) {
        returnType(double(0))
        if (x == 0) {
            # Probability of zero from both components
            # p_nb(0) = (r / (r + mu))^r
            p_nb_zero <- pow(r / (r + mu), r)
            prob <- psi + (1 - psi) * p_nb_zero
        } else {
            prob <- (1 - psi) * dnbinom(x, size = r, mu = mu, log = 0)
        }
        if (log) return(log(max(1.0e-30, prob))) else return(prob)
    }
)

#' @export
rzinb_because <- nimble::nimbleFunction(
    run = function(
        n = integer(0),
        mu = double(0),
        r = double(0),
        psi = double(0)
    ) {
        returnType(double(0))
        if (runif(1) < psi) {
            return(0)
        } else {
            return(rnbinom(1, size = r, mu = mu))
        }
    }
)

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
    # Vague prior is safe for Gaussian
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
        "a module is required (e.g., because.occupancy)."
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
    tau_var <- paste0("tau_e_", response, suffix)

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
            dzip_because = dzip_because,
            rzip_because = rzip_because
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
            dzinb_because = dzinb_because,
            rzinb_because = rzinb_because
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
            dnb_because = dnb_because,
            rnb_because = rnb_because
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
