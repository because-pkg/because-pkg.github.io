#' Run a NumPyro model via because_py
#'
#' @param eq_strings Character vector of equations
#' @param flat_data Flattened list of data
#' @param family Named character vector of families
#' @param priors Named list of priors
#' @param py_structures Processed structural matrices/functions for NumPyro
#' @param n_chains Number of MCMC chains
#' @param n_iter Number of iterations
#' @param n_warmup Number of warmup iterations
#' @param adapt_delta Target acceptance probability
#' @param max_treedepth Maximum tree depth
#' @param prior_scale_fixed Scale factor for fixed effects
#' @param quiet Logical, suppress output
#' @return Raw result from because_py$fit_numpyro_model
#' @export
run_numpyro_model <- function(eq_strings, flat_data, family, priors, py_structures,
                              n_chains, n_iter, n_warmup, adapt_delta, max_treedepth,
                              prior_scale_fixed, quiet) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required to use engine = 'numpyro'")
  }
  
  because_py <- tryCatch({
    reticulate::import("because.api")
  }, error = function(e) {
    stop("Failed to import 'because.api'. Make sure because_py is installed and reticulate is configured.")
  })
  
  py_result <- because_py$fit_numpyro_model(
    equations = eq_strings,
    data = flat_data,
    family = family,
    priors = priors,
    n_chains = as.integer(n_chains),
    n_samples = as.integer(n_iter),
    n_warmup = as.integer(n_warmup),
    cor_matrices = py_structures,
    adapt_delta = adapt_delta,
    max_treedepth = as.integer(max_treedepth),
    prior_scale_fixed = prior_scale_fixed
  )
  
  return(py_result)
}

#' Convert NumPyro output into coda mcmc.list
#'
#' @param py_result Result list from fit_numpyro_model
#' @return coda::mcmc.list
#' @export
format_numpyro_samples <- function(py_result) {
  raw_samples <- py_result$samples
  
  if (length(raw_samples) == 0) return(NULL)
  
  # Determine chains and iterations from the first parameter's dimensions
  first_param <- raw_samples[[1]]
  num_chains <- dim(first_param)[1]
  num_iters <- dim(first_param)[2]
  
  chain_list <- list()
  for (ch in 1:num_chains) {
    chain_mat <- NULL
    for (param_name in names(raw_samples)) {
      param_data <- raw_samples[[param_name]]
      
      if (length(dim(param_data)) == 2) {
        # Scalar parameter: shape (chains, iters)
        col <- matrix(param_data[ch, ], ncol = 1)
        colnames(col) <- param_name
        chain_mat <- if (is.null(chain_mat)) col else cbind(chain_mat, col)
      } else if (length(dim(param_data)) == 3) {
        # Vector parameter: shape (chains, iters, length)
        param_len <- dim(param_data)[3]
        cols <- matrix(param_data[ch, , ], ncol = param_len)
        colnames(cols) <- paste0(param_name, "[", 1:param_len, "]")
        chain_mat <- if (is.null(chain_mat)) cols else cbind(chain_mat, cols)
      }
    }
    chain_list[[ch]] <- coda::mcmc(chain_mat)
  }
  
  return(coda::mcmc.list(chain_list))
}
