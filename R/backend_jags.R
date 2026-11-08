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
  ens_ntree <- data[[".ensemble_ntree"]]
  ens_cpt   <- data[[".ensemble_n_chains"]] %||% n.chains %||% 2L
  ens_k_var <- data[[".ensemble_k_var"]] %||% "K_tree_multiPhylo"
  data[[".ensemble_ntree"]] <- NULL
  data[[".ensemble_n_chains"]] <- NULL
  data[[".ensemble_k_var"]] <- NULL

  if (!is.null(ens_ntree) && ens_ntree > 1L) {
    if (!is.null(extension_inits)) {
      extension_inits[[ens_k_var]] <- NULL
    }
    total_ens_chains <- as.integer(ens_ntree * ens_cpt)
    base_seed <- sample.int(1e6L, 1L)

    if (parallel && n.cores > 1L && total_ens_chains > 1L) {
      n_workers <- as.integer(min(n.cores, total_ens_chains))
      if (!quiet) {
        message(sprintf(
          "  [MultiPhylo Ensemble] Running %d JAGS chains (%d trees x %d chains/tree) in parallel across %d CPU cores...",
          total_ens_chains, ens_ntree, ens_cpt, n_workers
        ))
      }
      if (is.null(cl)) {
        cl <- parallel::makeCluster(n_workers)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      }
      run_ens_jags_chain <- function(idx, model_file, data, ens_k_var, ens_cpt, extension_inits, base_seed, monitor, n.burnin, n.iter, n.thin, n.adapt, quiet, model_string) {
        loadNamespace("rjags")
        k <- as.integer(ceiling(idx / ens_cpt))
        data_k <- data
        data_k[[ens_k_var]] <- k
        inits_one <- list(c(
          extension_inits,
          list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = base_seed + idx)
        ))
        model_k <- run_jags_model(model_file, data_k, inits_one, 1L, n.adapt, TRUE, model_string)
        if (n.burnin > 0) update(model_k, n.iter = n.burnin)
        samples_k <- sample_jags_model(model_k, monitor, n.iter, n.burnin, n.thin)
        mat_c <- as.matrix(samples_k[[1]])
        if (!ens_k_var %in% colnames(mat_c)) {
          k_mat <- matrix(as.numeric(k), nrow = nrow(mat_c), ncol = 1, dimnames = list(NULL, ens_k_var))
          mat_c <- cbind(mat_c, k_mat)
        }
        list(mcmc = coda::mcmc(mat_c), model = model_k)
      }
      parallel::clusterExport(cl, c("run_ens_jags_chain", "run_jags_model", "sample_jags_model"), envir = environment())
      parallel::clusterEvalQ(cl, {
        if (requireNamespace("because", quietly = TRUE)) library(because)
      })
      par_res <- parallel::parLapply(
        cl,
        seq_len(total_ens_chains),
        run_ens_jags_chain,
        model_file = model_file,
        data = data,
        ens_k_var = ens_k_var,
        ens_cpt = ens_cpt,
        extension_inits = extension_inits,
        base_seed = base_seed,
        monitor = monitor,
        n.burnin = n.burnin,
        n.iter = n.iter,
        n.thin = n.thin,
        n.adapt = n.adapt,
        quiet = quiet,
        model_string = model_string
      )
      all_chain_samples <- lapply(par_res, `[[`, "mcmc")
      last_model <- par_res[[length(par_res)]]$model
    } else {
      all_chain_samples <- list()
      last_model <- NULL
      for (k in seq_len(ens_ntree)) {
        if (!quiet) {
          message(sprintf("  [MultiPhylo Ensemble: Tree %d/%d] Sampling %d JAGS chain(s)...", k, ens_ntree, ens_cpt))
        }
        data_k <- data
        data_k[[ens_k_var]] <- as.integer(k)
        inits_k <- lapply(seq_len(ens_cpt), function(c_idx) {
          c(
            extension_inits,
            list(
              .RNG.name = "base::Wichmann-Hill",
              .RNG.seed = base_seed + (k - 1L) * ens_cpt + c_idx
            )
          )
        })
        model_k <- run_jags_model(model_file, data_k, inits_k, ens_cpt, n.adapt, quiet, model_string)
        if (n.burnin > 0) {
          update(model_k, n.iter = n.burnin)
        }
        samples_k <- sample_jags_model(model_k, monitor, n.iter, n.burnin, n.thin)
        for (c_idx in seq_len(ens_cpt)) {
          mat_c <- as.matrix(samples_k[[c_idx]])
          if (!ens_k_var %in% colnames(mat_c)) {
            k_mat <- matrix(as.numeric(k), nrow = nrow(mat_c), ncol = 1, dimnames = list(NULL, ens_k_var))
            mat_c <- cbind(mat_c, k_mat)
          }
          all_chain_samples[[length(all_chain_samples) + 1L]] <- coda::mcmc(mat_c)
        }
        last_model <- model_k
      }
    }

    return(list(
      samples = coda::mcmc.list(all_chain_samples),
      model = last_model,
      DIC = FALSE,
      WAIC = WAIC
    ))
  }

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

    base_seed <- sample.int(1e6L, 1L)

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
          .RNG.seed = base_seed + chain_id
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

    parallel::clusterExport(cl, c("run_single_chain", "run_jags_model", "sample_jags_model", "base_seed"), envir = environment())
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

    base_seed <- sample.int(1e6L, 1L)
    inits_list <- lapply(1:n.chains, function(i) {
      c(
        extension_inits,
        list(
          .RNG.name = "base::Wichmann-Hill",
          .RNG.seed = base_seed + i
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
