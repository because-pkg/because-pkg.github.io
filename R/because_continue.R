#' Continue MCMC Sampling for a Fitted Model
#'
#' Extends an existing \code{because} model by running additional MCMC
#' iterations, then (optionally) combining the new samples with the original
#' ones.  This is useful when traceplots or Rhat values suggest that chains
#' have not converged and more iterations are needed — without discarding the
#' burn-in and samples already obtained.
#'
#' @param object A fitted \code{because} model (class \code{"because"}).
#'   Must have been run with \code{dsep = FALSE} (i.e. a full MCMC fit, not
#'   a d-separation-only run).
#' @param n.iter Integer. Number of \emph{additional} post-warmup iterations
#'   to draw per chain, in the same units as the \code{n.iter} argument of
#'   \code{because()}.  The number of stored samples added per chain is
#'   \code{n.iter \%/\% thin}, where \code{thin} is inherited from the
#'   original run.
#' @param combine Logical (default \code{TRUE}). If \code{TRUE} the new
#'   samples are appended to the existing ones so the returned object contains
#'   all iterations.  If \code{FALSE} only the new samples are kept.
#' @param parallel Logical or \code{NULL} (default). \code{NULL} reuses the
#'   same setting as the original \code{because()} call (stored in
#'   \code{object$parallel}).  Set explicitly to \code{TRUE} or \code{FALSE}
#'   to override.  Only relevant for JAGS rebuild paths; ignored for Nimble
#'   (see Details).
#' @param n.cores Integer or \code{NULL}. Number of cores for parallel
#'   execution. \code{NULL} reuses the value from the original run.
#' @param n.adapt Integer. Adaptation iterations used when the JAGS model
#'   must be rebuilt from scratch (ignored when the live model object is
#'   reused, and not applicable to Nimble). Default \code{200}.
#' @param quiet Logical. Suppress progress messages. Default \code{FALSE}.
#'
#' @return The updated \code{because} object with extended (or replaced)
#'   samples and a recalculated summary including Rhat.  The live model
#'   object (\code{result$model} for JAGS, \code{result$nimble_compiled} for
#'   Nimble) is updated so that \code{because_continue()} can be called again
#'   without recompilation.
#'
#' @details
#' \strong{Engine support and practical guidance}
#'
#' \code{because_continue()} is most useful with the \strong{JAGS} engine,
#' which is the slowest of the three supported engines.  The function
#' provides a \emph{warm-start rebuild} path that works even across R session
#' restarts, making it a practical tool for extending long-running JAGS fits
#' without discarding the existing samples or redoing burn-in.
#'
#' \tabular{lll}{
#'   \strong{Engine} \tab \strong{When it works} \tab \strong{What happens} \cr
#'   \code{jags} \tab Always (see strategies below) \tab No or minimal burn-in \cr
#'   \code{nimble} \tab Same R session + \code{parallel = FALSE} only \tab \code{run(reset=FALSE)}, no recompilation \cr
#'   \code{numpyro} \tab Not supported \tab Re-run \code{because()} with larger \code{n.iter} \cr
#' }
#'
#' \strong{JAGS — sampling strategies (tried in order)}
#' \enumerate{
#'   \item \strong{Live model (sequential runs, same session)} — if
#'     \code{object$model} is still a valid \code{jags} object,
#'     \code{rjags::coda.samples()} is called directly.  No recompilation,
#'     no burn-in, exact continuation of all chains.
#'   \item \strong{Parallel rebuild} — if \code{parallel = TRUE}, each chain
#'     is recompiled in its own worker process and warm-started from the last
#'     sample of that chain (\code{alpha_*} and \code{beta_*} parameters
#'     only, since derived nodes such as \code{sigma_*_res} or
#'     \code{lambda_*} cannot be initialised directly in JAGS).
#'   \item \strong{Sequential rebuild} — a single \code{jags.model()} call
#'     with warm-start inits.  Automatically falls back to progressively
#'     safer init sets (alpha/beta only, then cold start) when deterministic
#'     nodes are encountered.
#' }
#'
#' The rebuild paths work \emph{across R sessions} because the JAGS model
#' code is stored in \code{object$model_code} and recompilation is fast
#' (seconds, not minutes).
#'
#' \strong{Nimble — live compiled object path}
#'
#' For Nimble models compiled in the \emph{current} R session with
#' \code{parallel = FALSE}, the C++ compiled MCMC object
#' (\code{object$nimble_compiled}) is reused directly via
#' \code{compiled_mcmc$run(reset = FALSE)}.  Each chain is reinitialised to
#' the last sample of that chain in the stored \code{mcmc.list}, so the
#' continuation starts from the exact posterior position — no recompilation
#' (which can take several minutes for large models).
#'
#' Nimble continuation is \strong{not available} when:
#' \itemize{
#'   \item The original run used \code{parallel = TRUE} — compiled objects
#'     live in worker processes and are destroyed when sampling finishes.
#'   \item R has been restarted — C++ compiled objects are session-specific
#'     and cannot be saved to disk.
#' }
#' In those cases, re-run \code{because()} with a larger \code{n.iter}.
#'
#' \strong{NumPyro}
#'
#' NumPyro uses HMC/NUTS, which converges so efficiently that extending an
#' existing run is rarely necessary.  Re-run \code{because()} with a larger
#' \code{n.iter} if more samples are needed.
#'
#' @examples
#' \dontrun{
#' # --- JAGS: sequential run, same session ---
#' # Live model object is reused — instant continuation, no recompilation.
#' fit <- because(equations = list(y ~ x), data = dat,
#'                n.iter = 2000, n.burnin = 500, n.chains = 3)
#' summary(fit)  # Rhat looks borderline
#' fit <- because_continue(fit, n.iter = 5000)
#' summary(fit)  # now 6500 effective iterations
#'
#' # --- JAGS: after restarting R ---
#' # Model is rebuilt from object$model_code with warm-start inits.
#' # parallel = NULL reuses the original setting automatically.
#' load("fit.RData")
#' fit <- because_continue(fit, n.iter = 5000)
#'
#' # --- Nimble: sequential run, same session ---
#' fit_nim <- because(equations = list(y ~ x), data = dat,
#'                    engine = "nimble", n.iter = 2000, n.burnin = 500,
#'                    parallel = FALSE)
#' fit_nim <- because_continue(fit_nim, n.iter = 5000)
#' }
#'
#' @export
because_continue <- function(
    object,
    n.iter   = 5000,
    combine  = TRUE,
    parallel = NULL,
    n.cores  = NULL,
    n.adapt  = 200,
    quiet    = FALSE
) {
  # ---- Validation ----
  if (!inherits(object, "because")) {
    stop("'object' must be a fitted 'because' model.")
  }
  if (is.null(object$samples)) {
    stop(paste(
      "'object' has no MCMC samples.",
      "Only full model fits (dsep = FALSE) can be continued."
    ))
  }

  # ---- Engine check ----
  engine <- if (!is.null(object$engine)) object$engine else "jags"

  if (engine == "numpyro") {
    stop(paste(
      "because_continue() does not yet support engine = 'numpyro'.",
      "NumPyro's HMC/NUTS is fast enough that re-running because() with a",
      "larger n.iter is the recommended approach."
    ))
  }

  # ---- Resolve parallel / n.cores ----
  if (is.null(parallel)) {
    parallel <- isTRUE(object$parallel)
  }
  if (is.null(n.cores)) {
    n.cores <- if (!is.null(object$n.cores)) object$n.cores
               else max(1L, parallel::detectCores() - 1L)
  }

  n_chains  <- length(object$samples)
  var_names <- colnames(as.matrix(object$samples[[1]]))
  thin_val  <- coda::thin(object$samples)

  # ---- Strategy N: Nimble continuation ----
  # Nimble compiled objects live in the current R session. We reset the
  # compiled model to the last sample of each chain in turn, then call
  # compiled_mcmc$run(reset = FALSE) — no recompilation needed.
  if (engine == "nimble") {
    if (is.null(object$nimble_compiled) || is.null(object$nimble_cmodel)) {
      stop(paste(
        "Live Nimble compiled objects not found (object$nimble_compiled /",
        "object$nimble_cmodel are NULL).",
        "This happens when the model was run with parallel = TRUE or in a",
        "previous R session. Please re-run because() with a larger n.iter."
      ))
    }

    if (!quiet) message(sprintf(
      "Continuing Nimble model via run(reset=FALSE) (%d chain(s), %d additional iterations)...",
      n_chains, n.iter
    ))

    new_chain_list <- vector("list", n_chains)

    for (ch in seq_len(n_chains)) {
      # --- reset compiled model to last sample of this chain ---
      last_row <- tail(as.matrix(object$samples[[ch]]), 1L)[1L, ]

      # Set scalar parameters
      scalar_params <- names(last_row)[!grepl("[", names(last_row), fixed = TRUE)]
      for (pname in scalar_params) {
        tryCatch(
          object$nimble_cmodel[[pname]] <- last_row[[pname]],
          error = function(e) {}   # skip read-only / deterministic nodes
        )
      }

      # Set array parameters (group columns by base name)
      array_cols <- names(last_row)[grepl("[", names(last_row), fixed = TRUE)]
      if (length(array_cols) > 0L) {
        base_names <- unique(sub("\\[.*\\]$", "", array_cols))
        for (bname in base_names) {
          cols  <- array_cols[startsWith(array_cols, paste0(bname, "["))]
          vals  <- last_row[cols]
          # parse indices and sort
          idx   <- as.integer(regmatches(cols,
                     regexpr("(?<=\\[)\\d+(?=\\])", cols, perl = TRUE)))
          vals  <- vals[order(idx)]
          tryCatch(
            object$nimble_cmodel[[bname]] <- unname(vals),
            error = function(e) {}
          )
        }
      }

      # Propagate to deterministic nodes
      tryCatch(object$nimble_cmodel$calculate(), error = function(e) {})

      # Run MCMC from current state, no reset
      object$nimble_compiled$run(
        niter    = n.iter,
        reset    = FALSE,
        resetMV  = TRUE    # reset sample store so mvSamples holds THIS run only
      )

      ch_mat <- as.matrix(object$nimble_compiled$mvSamples)
      # Keep only monitored parameters (exclude log_lik etc.)
      keep   <- colnames(ch_mat)[!grepl("^log_lik", colnames(ch_mat))]
      new_chain_list[[ch]] <- coda::mcmc(
        ch_mat[, keep, drop = FALSE],
        thin = thin_val
      )
    }

    new_samples <- coda::mcmc.list(new_chain_list)
    new_model   <- object$model   # R model object (unchanged)

    # Fall through to combine/summary block below
  }

  # ---- Strategy 1: reuse live rjags model object ----
  # Only applies to JAGS runs; skipped for Nimble (handled above).
  if (engine == "jags") {
  new_samples <- NULL
  new_model   <- NULL

  if (!is.null(object$model) && inherits(object$model, "jags")) {
    if (!quiet) message(sprintf(
      "Live JAGS model found. Sampling %d additional iterations across %d chain(s)...",
      n.iter, n_chains
    ))
    tryCatch({
      new_samples <- rjags::coda.samples(
        model          = object$model,
        variable.names = var_names,
        n.iter         = n.iter,
        thin           = thin_val
      )
      new_model <- object$model
    }, error = function(e) {
      if (!quiet) message(
        "Live model object unusable. Falling back to rebuild..."
      )
    })
  }

  # ---- Strategy 2 & 3: rebuild ----
  if (is.null(new_samples)) {
    model_code <- object$model_code
    if (is.null(model_code) || !nzchar(model_code)) {
      stop(paste(
        "No live model object and no model code found in object$model_code.",
        "Cannot continue this run."
      ))
    }
    if (is.null(object$data)) {
      stop("object$data is missing. Cannot rebuild the model.")
    }

    # Build warm-start init sets.
    # JAGS cannot initialise *deterministic* nodes (sigma_*_res, lambda_*,
    # sigma_total_*, etc.).  We therefore provide two levels:
    #   (a) alpha_* and beta_* — always stochastic, safe to initialise
    #   (b) empty list — cold start fallback
    safe_scalar <- var_names[
      grepl("^(alpha|beta)_", var_names) &
      !grepl("[", var_names, fixed = TRUE)
    ]

    make_chain_inits <- function(chain_idx, keep_vars) {
      row <- tail(as.matrix(object$samples[[chain_idx]]), 1L)[1L, ]
      as.list(row[names(row) %in% keep_vars])
    }

    model_file <- tempfile(fileext = ".jg")
    writeLines(model_code, model_file)
    on.exit(unlink(model_file), add = TRUE)

    # ---- Strategy 2: parallel rebuild ----
    if (parallel && n.cores > 1L && n_chains > 1L) {
      if (!quiet) message(sprintf(
        "Rebuilding %d chain(s) in parallel on %d core(s) (%d additional iterations each)...",
        n_chains, min(n.cores, n_chains), n.iter
      ))

      run_chain <- function(chain_id, model_file, data, inits,
                            var_names, n.iter, thin_val, n.adapt, quiet) {
        loadNamespace("rjags")
        m <- tryCatch(
          rjags::jags.model(
            file     = model_file,
            data     = data,
            inits    = c(inits, list(.RNG.name  = "base::Wichmann-Hill",
                                     .RNG.seed  = 12345L + chain_id)),
            n.chains = 1L,
            n.adapt  = n.adapt,
            quiet    = quiet
          ),
          error = function(e) {
            if (grepl("non-variable node|Cannot set value", e$message,
                      ignore.case = TRUE)) {
              rjags::jags.model(
                file     = model_file,
                data     = data,
                inits    = list(.RNG.name = "base::Wichmann-Hill",
                                .RNG.seed = 12345L + chain_id),
                n.chains = 1L,
                n.adapt  = n.adapt,
                quiet    = quiet
              )
            } else stop(e)
          }
        )
        s <- rjags::coda.samples(m, variable.names = var_names,
                                 n.iter = n.iter, thin = thin_val)
        list(samples = s, model = m)
      }

      cl <- parallel::makeCluster(min(n.cores, n_chains))
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, "run_chain", envir = environment())

      chain_inits <- lapply(seq_len(n_chains), make_chain_inits,
                            keep_vars = safe_scalar)

      chain_results <- parallel::parLapply(
        cl,
        seq_len(n_chains),
        function(i) {
          run_chain(i, model_file, data, chain_inits[[i]],
                    var_names, n.iter, thin_val, n.adapt, quiet)
        }
      )

      new_samples <- coda::mcmc.list(lapply(chain_results,
                                            function(x) x$samples[[1]]))
      new_model   <- chain_results[[1]]$model   # first chain, for future use

    } else {
      # ---- Strategy 3: sequential rebuild ----
      if (!quiet) message(sprintf(
        "Rebuilding JAGS model sequentially (%d chain(s), %d additional iterations)...",
        n_chains, n.iter
      ))

      all_inits  <- lapply(seq_len(n_chains), make_chain_inits,
                           keep_vars = safe_scalar)

      try_compile <- function(inits, label) {
        tryCatch(
          rjags::jags.model(
            file     = model_file,
            data     = object$data,
            inits    = inits,
            n.chains = n_chains,
            n.adapt  = n.adapt,
            quiet    = quiet
          ),
          error = function(e) {
            if (grepl("non-variable node|Cannot set value", e$message,
                      ignore.case = TRUE)) {
              if (!quiet) message(
                "  Warm-start (", label, ") hit a deterministic node. ",
                "Trying safer fallback..."
              )
              NULL
            } else stop(e)
          }
        )
      }

      new_model <- try_compile(all_inits, "alpha/beta warm-start")
      if (is.null(new_model)) {
        if (!quiet) message("  Falling back to cold start (no initial values).")
        new_model <- rjags::jags.model(
          file     = model_file,
          data     = object$data,
          n.chains = n_chains,
          n.adapt  = n.adapt,
          quiet    = quiet
        )
      }

      new_samples <- rjags::coda.samples(
        model          = new_model,
        variable.names = var_names,
        n.iter         = n.iter,
        thin           = thin_val
      )
    }
  }
  }  # end if (engine == "jags")

  # ---- Combine or replace ----
  if (combine) {
    if (!quiet) message("Combining new samples with existing samples...")
    combined <- coda::mcmc.list(
      mapply(
        function(old_ch, new_ch) {
          old_mat <- as.matrix(old_ch)
          new_mat <- as.matrix(new_ch)
          shared  <- intersect(colnames(old_mat), colnames(new_mat))
          coda::mcmc(
            rbind(old_mat[, shared, drop = FALSE],
                  new_mat[, shared, drop = FALSE]),
            thin = thin_val
          )
        },
        object$samples, new_samples,
        SIMPLIFY = FALSE
      )
    )
  } else {
    combined <- new_samples
  }

  # ---- Update result ----
  result         <- object
  result$samples <- combined
  result$model   <- new_model   # retain live JAGS object for future continuations
  # For Nimble: compiled objects are reference-class objects (R5), so they are
  # already updated in-place; we just make sure the slots are not dropped.
  if (engine == "nimble") {
    result$nimble_compiled <- object$nimble_compiled
    result$nimble_cmodel   <- object$nimble_cmodel
  }

  # Recalculate summary + Rhat
  sum_stats <- summary(combined)
  if (n_chains > 1L) {
    tryCatch({
      psrf     <- coda::gelman.diag(combined, multivariate = FALSE)$psrf
      common   <- intersect(rownames(sum_stats$statistics), rownames(psrf))
      rhat_vec <- rep(NA_real_, nrow(sum_stats$statistics))
      names(rhat_vec) <- rownames(sum_stats$statistics)
      rhat_vec[common] <- psrf[common, "Point est."]
      sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = rhat_vec)
    }, error = function(e) {
      if (!quiet) message("Note: Rhat calculation failed: ", e$message)
    })
  }
  result$summary <- sum_stats

  if (!quiet) {
    n_total  <- nrow(as.matrix(combined[[1L]]))
    n_added  <- n.iter %/% thin_val
    n_before <- n_total - n_added
    message(sprintf(
      "Done. Samples per chain: %d  (was %d, added %d)",
      n_total, n_before, n_added
    ))
  }

  return(result)
}
