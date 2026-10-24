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

#' Execute full JAGS pipeline (sequential or parallel) for because()
#'
#' @keywords internal
run_jags_pipeline <- function(
  model_file, model_string, data, extension_inits, monitor,
  n.chains, n.iter, n.burnin, n.thin, n.adapt,
  quiet, verbose, parallel, n.cores, cl = NULL,
  DIC = FALSE, WAIC = FALSE
) {
  if (parallel && n.cores > 1 && n.chains > 1) {
    # Parallel execution
    message(sprintf(
      "Running %d chains in parallel on %d cores...",
      n.chains,
      n.cores
    ))

    # Setup cluster if not provided
    if (is.null(cl)) {
      cl <- parallel::makeCluster(n.cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }

    # Helper function to run a single chain
    run_single_chain <- function(
      chain_id,
      model_file,
      data,
      monitor,
      n.burnin,
      n.iter,
      n.thin,
      n.adapt,
      quiet
    ) {
      if (!requireNamespace("rjags", quietly = TRUE)) {
        stop("Package 'rjags' is required for parallel execution.")
      }
      loadNamespace("rjags")

      inits_list <- c(
        extension_inits,
        list(
          .RNG.name = "base::Wichmann-Hill",
          .RNG.seed = 12345 + chain_id
        )
      )

      par_inits <- list(inits_list)
      model <- run_jags_model(model_file, data, par_inits, 1L, n.adapt, quiet, model_string)

      if (n.burnin > 0) {
        update(model, n.iter = n.burnin)
      }

      samples <- sample_jags_model(model, monitor, n.iter, n.burnin, n.thin)

      return(list(samples = samples, model = model))
    }

    parallel::clusterExport(cl, c("run_single_chain", "run_jags_model", "sample_jags_model"), envir = environment())
    parallel::clusterEvalQ(cl, {
      if (requireNamespace("because", quietly = TRUE)) library(because)
    })

    if (!quiet) {
      message(sprintf("Sampling %d chains in parallel...", n.chains))
    }

    chain_results <- parallel::parLapply(cl, seq_len(n.chains), function(i) {
      res <- run_single_chain(
        i,
        model_file,
        data,
        monitor,
        n.burnin,
        n.iter,
        n.thin,
        n.adapt,
        quiet
      )
      return(res)
    })

    if (!quiet) {
      message("All chains completed.")
    }

    if (!is.null(chain_results[[1]]$samples)) {
      samples <- coda::mcmc.list(lapply(chain_results, function(x) {
        x$samples[[1]]
      }))
    } else {
      samples <- NULL
    }

    model <- chain_results[[1]]$model
  } else {
    # Sequential execution (default)
    if (verbose) {
      message("--- JAGS MODEL STRING ---")
      message(model_string)
    }
    if (verbose) {
      cat(
        "\n--- DATA LIST NAMES ---\n",
        paste(names(data), collapse = ", "),
        "\n"
      )
    }

    inits_list <- lapply(1:n.chains, function(i) {
      c(
        extension_inits,
        list(
          .RNG.name = "base::Wichmann-Hill",
          .RNG.seed = 12345 + i
        )
      )
    })

    model <- run_jags_model(model_file, data, inits_list, n.chains, n.adapt, quiet, model_string)
    if (n.burnin > 0) {
      update(model, n.iter = n.burnin)
    }

    if (n.chains < 2 && (DIC || WAIC)) {
      warning(
        "DIC and WAIC require at least 2 chains. Disabling calculation."
      )
      DIC <- FALSE
      WAIC <- FALSE
    }

    samples <- sample_jags_model(model, monitor, n.iter, n.burnin, n.thin)
  }

  return(list(
    samples = samples,
    model = model,
    DIC = DIC,
    WAIC = WAIC
  ))
}
