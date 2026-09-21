#' Run a JAGS model with standard because settings
#'
#' @param model_file Path to the JAGS model file or text connection
#' @param data List of data for the model
#' @param inits_list List of initial values
#' @param n.chains Number of MCMC chains
#' @param n.adapt Number of adaptation steps
#' @param quiet Logical, suppress output
#' @param model_string The model string (for error reporting)
#' @return A compiled jags.model object
#' @export
run_jags_model <- function(model_file, data, inits_list, n.chains, n.adapt, quiet, model_string) {
  tryCatch(
    {
      rjags::jags.model(
        model_file,
        data = data,
        inits = inits_list,
        n.chains = n.chains,
        n.adapt = n.adapt,
        quiet = quiet
      )
    },
    error = function(e) {
      if (!quiet) {
        message("\nCRITICAL JAGS ERROR during compilation:")
        message(e$message)
        message("Check your model code syntax or data dimensions.\n")
      }
      stop(paste(e, "\n\n", model_string))
    }
  )
}

#' Sample from a compiled JAGS model
#'
#' @param model Compiled jags.model object
#' @param monitor Character vector of variables to monitor
#' @param n.iter Number of iterations
#' @param n.burnin Number of burnin iterations
#' @param n.thin Thinning interval
#' @return mcmc.list of samples
#' @export
sample_jags_model <- function(model, monitor, n.iter, n.burnin, n.thin) {
  if (n.iter > n.burnin) {
    rjags::coda.samples(
      model,
      variable.names = monitor,
      n.iter = n.iter - n.burnin,
      thin = n.thin
    )
  } else {
    stop("n.iter must be greater than n.burnin")
  }
}
