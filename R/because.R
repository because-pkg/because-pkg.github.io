

#' @title Run a Bayesian Structural Equation Model (Because)
#'
#' @description
#' Fits a Bayesian Structural Equation Model (SEM) using JAGS, NIMBLE, or NumPyro.
#' Supports multiscale (hierarchical) data, custom covariance structures
#' (phylogenetic, spatial, etc.), missing data imputation, and d-separation
#' global fit testing.
#'
#' @param equations A list of R formulas describing the structural model. Each formula represents 
#'    a causal path (e.g., `Y ~ X1 + X2`). Categorical predictors are automatically handled.
#' @param data A data frame containing all variables. For **multiscale models**, 
#'   this must be a named list of data frames (one per scale/level).
#' @param family A named character vector specifying the distribution for each response 
#'   (e.g., `c(Y = "poisson", M = "gaussian")`). Supports "gaussian" (default), "poisson", 
#'   "binomial", "bernoulli", "ordinal", "multinomial", "negbinomial", "zip", and "zinb".
#' @param dsep Logical; if `TRUE`, performs d-separation tests to evaluate the global 
#'   fit of the DAG against the data. Highly recommended for causal validation.
#' @param dsep_max_obs (Diagnostic Optimization) Maximum number of observations used 
#'   for any single d-separation test (default = 10,000). For massive datasets, 
#'   this provides a significant speedup without losing structural diagnostic 
#'   power. Set to \code{Inf} to disable subsampling.
#' @param random Optional formula for **global random intercepts** (e.g., `~(1|Site)`). 
#'   Deprecated in favor of inline specification: `Y ~ X + (1|Site)`.
#' @param latent Optional character vector of latent (unmeasured) variables.
#' @param id_col Character string for the column identifying units (e.g. species names).
#' @param structure Optional covariance structure (e.g., a phylogenetic tree or spatial matrix) 
#'   for modeling correlated residuals.
#' @param multiscale (Multiscale/Hierarchical) A string defining the nesting structure 
#'   (e.g., `"site > individual"`). Semicolons separate parallel branches 
#'   (e.g., `"site > obs; species > obs"`). Formerly called `hierarchy`.
#' @param levels (Hierarchical) Optional named list mapping variables to their home levels. 
#'   If omitted, `because` will attempt to auto-detect levels based on column availability.
#' @param link_vars (Hierarchical) A named character vector specifying the columns used to 
#'   link data frames across levels (e.g., `c(site = "SiteID")`).
#' @param engine Bayesian backend: `"jags"` (default), `"nimble"`, or `"numpyro"`. 
#'   Note: For Gaussian models, the framework now uses a robust, scale-invariant `Uniform(0, 100)`
#'   prior on the standard deviation across all engines, replacing legacy precision priors.
#' @param n.iter Total MCMC iterations per chain (default = 12500).
#' @param n.burnin Number of burn-in iterations (default = 20% of `n.iter`).
#' @param n.chains Number of independent MCMC chains (default = 3).
#' @param parallel Logical; if `TRUE`, runs MCMC chains in parallel.
#' @param n.cores Number of CPU cores for parallel execution.
#' @param monitor Monitoring mode (e.g. "interpretable", "all", or "minimal").
#' @param nimble_samplers Optional named list of user-specified samplers for NIMBLE.
#' @param n.thin Thinning interval for MCMC chains (default = 10).
#' @param DIC Logical; if `TRUE`, calculates the Deviance Information Criterion (JAGS only).
#' @param WAIC Logical; if `TRUE`, calculates the Watanabe-Akaike Information Criterion.
#' @param n.adapt Number of iterations for the adaptation phase.
#' @param adapt_delta Target acceptance probability for the NUTS sampler (NumPyro only, default = 0.95).
#'   Increase towards 1 (e.g., `0.99`) for complex posteriors with funnel geometry or many
#'   competing variance components. Higher values slow sampling but improve mixing.
#' @param prior_scale_fixed Numeric; scale parameter for Cauchy priors on fixed effects. Defaults to `sqrt(2)/2`.
#' @param max_treedepth Maximum tree depth for the NUTS sampler (NumPyro only, default = 10).
#'   Increase to 12 or 14 if you see many divergent transitions or very low n.eff.
#' @param quiet Logical; if `TRUE`, suppresses status messages and progress bars.
#' @param verbose Logical; if `TRUE`, prints detailed debugging information.
#' @param variability Optional named character vector specifying the type of residual variance (e.g., "fixed", "random").
#' @param distribution Alias for `family` (deprecated).
#' @param latent_method Method for handling latent variables ("correlations" for MAG or "explicit").
#' @param standardize_latent Logical; if `TRUE`, standardizes latent variables to unit variance.
#' @param fix_latent How to fix the scale of latent variables ("loading" or "variance").
#' @param cl Optional pre-existing cluster for parallel execution.
#' @param ic_recompile Logical; if `TRUE`, recompiles the model for DIC/WAIC calculation.
#' @param hierarchy Alias for `multiscale` (deprecated).
#' @param fix_residual_variance Optional named vector fixing residual variances to known values.
#' @param priors Optional list of custom priors for model parameters.
#' @param reuse_models Optional list of previously fitted models to speed up d-separation tests.
#' @param expand_ordered Logical; if `TRUE`, expands ordered factors using monotonic effects.
#' @param structure_multi Logical; if `TRUE`, handles multiple covariance structures.
#' @param structure_levels Mapping of variables to covariance structure levels.
#' @param aggregate_crossscale Optional. Set to `"all"` to aggregate all cross-scale d-sep tests (evaluating coarse-scale conditional independencies by aggregating fine-scale data via group means), or provide a numeric vector specifying specific test indices to aggregate. By default, paths between a fine-scale and coarse-scale variable are flattened with hierarchical random effects instead.
#' @param ... Additional arguments passed to the underlying model engines.
#'
#' @return An object of class \code{"because"} containing:
#'   \item{samples}{MCMC samples (mcmc.list).}
#'   \item{parameter_map}{Data frame mapping parameter names to model variables.}
#'   \item{model}{JAGS model code.}
#'   \item{summary}{Summary of posterior samples.}
#'   \item{dsep_results}{List of d-separation test results (if dsep=TRUE).}
#'   \item{DIC}{Deviance Information Criterion (if DIC=TRUE).}
#'   \item{WAIC}{Watanabe-Akaike Information Criterion (if WAIC=TRUE).}
#'
#' @examples
#' \dontrun{
#' # 1. Simple Path Model
#' df <- data.frame(Y = rnorm(100), M = rnorm(100), X = rnorm(100))
#' eqs <- list(M ~ X, Y ~ M + X)
#' fit <- because(eqs, data = df, dsep = TRUE)
#' summary(fit)
#'
#' # 2. Multiscale Model (Auto-detection)
#' # Suppose we have year-level data and individual-level data
#' year_df <- data.frame(Year = 1:5, Temp = rnorm(5))
#' ind_df  <- data.frame(Year = rep(1:5, each=10), Mass = rnorm(50))
#' data_list <- list(yr = year_df, ind = ind_df)
#' 
#' # Mass depends on Temp (cross-scale link via Year)
#' fit_h <- because(
#'   equations  = list(Mass ~ Temp),
#'   data       = data_list,
#'   multiscale = "yr > ind",
#'   link_vars  = c(yr = "Year")
#' )
#'
#' # 3. Random Intercepts (lme4-style)
#' # Note: Grouping variables (Site) must be in the data
#' df$Site <- rep(1:10, each=10)
#' fit_re <- because(list(Y ~ X + (1|Site)), data = df)
#' }
#'
#' @export
#' @import coda
#' @import methods
#' @importFrom rjags jags.model coda.samples dic.samples jags.samples
#' @importFrom stats formula terms setNames start var na.omit update
#' @importFrom utils capture.output head tail data flush.console
#' @importFrom coda gelman.diag effectiveSize
because <- function(
  equations,
  data,
  id_col = NULL,
  structure = NULL,
  engine = "jags",
  monitor = "interpretable",
  nimble_samplers = NULL,
  n.chains = NULL,
  n.iter = NULL,
  n.burnin = NULL,
  n.thin = NULL,
  adapt_delta = 0.95,
  max_treedepth = 10,
  DIC = TRUE,
  WAIC = FALSE,
  n.adapt = NULL,
  quiet = FALSE,
  verbose = FALSE,
  dsep = FALSE,
  dsep_max_obs = 10000,
  variability = NULL,
  family = NULL,
  distribution = NULL,
  latent = NULL,
  latent_method = "correlations",
  standardize_latent = TRUE,
  fix_latent = "loading",
  prior_scale_fixed = NULL,
  parallel = FALSE,
  n.cores = parallel::detectCores() - 1,
  cl = NULL,
  ic_recompile = FALSE,
  random = NULL,
  levels = NULL,
  multiscale = NULL,
  hierarchy = NULL,
  link_vars = NULL,
  fix_residual_variance = NULL,
  priors = NULL,
  reuse_models = NULL,
  expand_ordered = FALSE,
  structure_multi = NULL,
  structure_levels = NULL,
  aggregate_crossscale = NULL,
  ...
) {
  original_call <- match.call(expand.dots = TRUE)
  args <- list(...)
  
  if (!is.null(link_vars)) {
      link_vars <- as.list(link_vars)
  }

  # Un-nest equations if passed as nested list e.g. list(list(...))
  if (is.list(equations) && length(equations) == 1 && is.list(equations[[1]]) && !inherits(equations[[1]], "formula")) {
    equations <- equations[[1]]
  }

  # --- Handle Terminology Aliases (Backward Compatibility) ---
  # hierarchy -> multiscale
  if (is.null(multiscale) && !is.null(hierarchy)) multiscale <- hierarchy
  if (is.null(hierarchy) && !is.null(multiscale)) hierarchy <- multiscale
  
  # --- Handle Deprecated Arguments ---
  # tree -> structure
  if (is.null(structure) && !is.null(args$tree)) {
    warning("Argument 'tree' is deprecated. Please use 'structure' instead.")
    structure <- args$tree
  }
  
  # multi.tree -> is_multi (passed to because_model)
  is_multi <- if (!is.null(args$is_multi)) args$is_multi else if (!is.null(args$multi.tree)) args$multi.tree else FALSE
  if (!is.null(args$multi.tree)) {
    warning("Argument 'multi.tree' is deprecated. Please use 'is_multi' instead.")
  }

  # Allow string input (e.g., multiline equations)
  if (is.character(equations) && length(equations) == 1) {
    # Split by newline or semicolon
    eq_lines <- unlist(strsplit(equations, "[\n;]"))
    eq_lines <- trimws(eq_lines)
    eq_lines <- eq_lines[eq_lines != ""]
    
    # Convert to formulas
    equations <- lapply(eq_lines, stats::as.formula)
  }

  # Allow single formula input
  if (inherits(equations, "formula")) {
    equations <- list(equations)
  }

  engine <- match.arg(tolower(engine), c("jags", "nimble", "numpyro"))
  
  # Auto-adjust MCMC defaults for NUTS sampler (NumPyro) if not explicitly provided
  if (engine == "numpyro") {
    if (is.null(n.iter))   n.iter   <- 2000
    if (is.null(n.thin))   n.thin   <- 1
    if (is.null(n.burnin)) n.burnin <- floor(n.iter / 2)

  }
  if (engine == "numpyro") {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("The 'reticulate' package is required when engine = 'numpyro'.")
    }
    # Automatically try to bind to the default r-reticulate virtual environment
    # before attempting to import the module, to save the user from doing it manually.
    tryCatch({
      # Bind to because_env with required = TRUE so we never silently fall
      # back to a different environment that may be missing because_py.
      if (reticulate::virtualenv_exists("because_env")) {
        reticulate::use_virtualenv("because_env", required = TRUE)
      } else if (reticulate::condaenv_exists("because_env")) {
        reticulate::use_condaenv("because_env", required = TRUE)
      }
    }, error = function(e) NULL)
    # ------------------------------------------------------------------
    # Thread-limit setup: must happen BEFORE reticulate::import("because.api")
    # because api.py now performs a lazy "import jax" on first fit() call,
    # and JAX reads thread counts from the C++ runtime at that moment.
    #
    # Strategy:
    #   1. Set vars in R's environment (covers child processes spawned by R).
    #   2. Mirror into Python's os.environ in case Python is already running
    #      (reticulate shares one Python process across the session).
    #   3. Set XLA_FLAGS for the correct device count.
    # ------------------------------------------------------------------
    target_cores <- as.integer(if (parallel) min(n.cores, n.chains) else 1L)

    # --- R-level env vars ---
    .thread_vars <- list(
      OMP_NUM_THREADS            = "1",
      OPENBLAS_NUM_THREADS       = "1",
      GOTO_NUM_THREADS           = "1",
      MKL_NUM_THREADS            = "1",
      MKL_DOMAIN_NUM_THREADS     = "1",
      NUMEXPR_NUM_THREADS        = "1",
      LLVM_NUM_THREADS           = "1",
      TF_NUM_INTEROP_THREADS     = as.character(target_cores),
      TF_NUM_INTRAOP_THREADS     = as.character(target_cores),
      XLA_PYTHON_CLIENT_PREALLOCATE = "false"
    )
    # Only set vars not already set by the user
    for (.v in names(.thread_vars)) {
      if (!nzchar(Sys.getenv(.v))) {
        do.call(Sys.setenv, stats::setNames(list(.thread_vars[[.v]]), .v))
      }
    }

    # XLA_FLAGS: merge device-count and Eigen flags without overwriting user flags
    .current_xla <- Sys.getenv("XLA_FLAGS")
    .xla_additions <- character(0)
    if (!grepl("--xla_force_host_platform_device_count", .current_xla))
      .xla_additions <- c(.xla_additions,
                          paste0("--xla_force_host_platform_device_count=", target_cores))
    if (!grepl("--xla_cpu_multi_thread_eigen", .current_xla))
      .xla_additions <- c(.xla_additions, "--xla_cpu_multi_thread_eigen=false")
    if (!grepl("intra_op_parallelism_threads", .current_xla))
      .xla_additions <- c(.xla_additions, paste0("intra_op_parallelism_threads=", target_cores))
    if (!grepl("inter_op_parallelism_threads", .current_xla))
      .xla_additions <- c(.xla_additions, paste0("inter_op_parallelism_threads=", target_cores))
    if (length(.xla_additions) > 0)
      Sys.setenv(XLA_FLAGS = trimws(paste(.current_xla, paste(.xla_additions, collapse = " "))))

    # --- Mirror vars into Python's os.environ (if Python is already running) ---
    if (reticulate::py_available(initialize = FALSE)) {
      .py_set_code <- paste(
        "import os",
        paste(sapply(names(.thread_vars), function(.v) {
          sprintf("os.environ.setdefault('%s', '%s')", .v, .thread_vars[[.v]])
        }), collapse = "\n"),
        sprintf("os.environ.setdefault('XLA_FLAGS', '%s')", Sys.getenv("XLA_FLAGS")),
        sep = "\n"
      )
      tryCatch(
        reticulate::py_run_string(.py_set_code),
        error = function(e) NULL
      )
    }

    tryCatch({
      because_py <- reticulate::import("because.api")
    }, error = function(e) {
      current_env <- "unknown"
      tryCatch({
        current_env <- reticulate::py_config()$python
      }, error = function(e) {})
      stop(sprintf(paste0(
        "Failed to import python module 'because.api'.\n",
        "Python is currently running from: %s\n\n",
        "This usually means Python was initialized to a different environment\n",
        "before 'library(because)' was called (e.g. by RStudio or another package).\n\n",
        "Quick fix --- add this line BEFORE library(because) in your script:\n",
        "  reticulate::use_virtualenv('because_env', required = TRUE)\n\n",
        "If because_env does not exist yet, install it first with:\n",
        "  install_because_numpyro()"
      ), current_env))
    })
  }

  # --- Robust Parameter Initialization (Handle NULL inputs from recursive calls) ---
  if (is.null(n.chains)) n.chains <- 3
  if (is.null(n.iter)) n.iter <- 12500
  if (is.null(n.thin)) n.thin <- 10
  if (is.null(n.burnin)) n.burnin <- floor(n.iter / 5)
  if (is.null(n.adapt)) n.adapt <- floor(n.iter / 5)

  # [FIX] Ensure all MCMC parameters are coercible to non-negative integer for JAGS
  n.chains <- as.integer(max(1, n.chains))
  n.iter <- as.integer(max(0, n.iter))
  n.thin <- as.integer(max(1, n.thin))
  n.burnin <- as.integer(max(0, n.burnin))
  n.adapt <- as.integer(max(0, n.adapt))

  if (!quiet) {
    is_list_debug <- is.list(data) && !is.data.frame(data)
  }

  # [PURGE] Legacy tree-to-structure alias removed. All structures handled via 'structure' argument.

  # [PURGE] Legacy structure check moved to top-level argument handler.

  # Validate inputs
  # Input validation
  if (is.null(data)) {
    stop("Argument 'data' must be provided.")
  }
  if (is.null(equations)) {
    stop("Argument 'equations' must be provided.")
  }

  # --- Clean Data & Family Normalization ---
  prep_res <- preprocess_because_args(data, equations, distribution, family, structure)
  data <- prep_res$data
  scale_info <- prep_res$scale_info
  family <- prep_res$family
  family_objects <- prep_res$family_objects
  family_obj <- prep_res$family_obj
  structure_obj <- prep_res$structure_obj
  equations <- prep_res$equations

  # --- Auto-Stacking (Multispecies Input) ---
  # If data is a list containing species-specific matrices (and equation is generic Y ~ ...),
  # convert to a single stacked dataframe (Long Format) with (1|SpeciesID).
  stack_res <- auto_stack_multispecies_data(data, equations, quiet = quiet)
  if (stack_res$is_stacked) {
    data <- stack_res$data

    # Append random effect for SpeciesID
    if (!is.null(stack_res$random_part)) {
      if (is.null(random)) {
        random <- as.formula(paste("~", stack_res$random_part))
      } else {
        r_str <- deparse(random)
        if (grepl("^~", r_str)) {
          r_str <- sub("^~", "", r_str)
        }
        # Avoid duplicating if user already added it?
        if (!grepl("SpeciesID", r_str)) {
          random <- as.formula(paste("~", r_str, stack_res$random_part))
        }
      }
      if (!quiet) {
        message("Auto-Stacking added random effect: ", stack_res$random_part)
      }
    }
  }

  # --- Selective Variable Collection ---
  # To avoid memory issues with large datasets, we only process variables needed by the model.

  # 1. Identify all variables used in equations (fixed + random)
  all_eq_vars <- unique(unlist(lapply(equations, all.vars)))

  # 2. Identify variables in global random effects
  global_random_vars <- if (!is.null(random)) {
    if (inherits(random, "formula")) {
      all.vars(random)
    } else {
      unique(unlist(lapply(random, all.vars)))
    }
  } else {
    character(0)
  }

  # 3. Identify categorical predictors that need dummy variables
  # These are variables on the RHS of fixed-effect formulas
  parsed_random_temp <- extract_random_effects(equations)
  fixed_eqs_temp <- parsed_random_temp$fixed_equations
  fixed_predictors <- unique(unlist(lapply(fixed_eqs_temp, function(eq) {
    if (length(eq) == 3) all.vars(eq[[3]]) else character(0)
  })))


  # 4. Combine all needed variables
  model_vars <- unique(c(
    all_eq_vars,
    global_random_vars,
    id_col,
    link_vars, # [HIERARCHY FIX] Keep link variables for hierarchical data prep
    if (!is.null(family)) names(family),
    if (!is.null(variability) && !is.character(variability)) names(variability),
    "N" # Explicitly keep N if provided
  ))

  # 5. Filter data frame early to save memory if it's a data.frame
  if (is.data.frame(data) || (is.list(data) && !is.data.frame(data))) {
    available_vars <- intersect(names(data), model_vars)
    if (length(available_vars) > 0) {

      # [NEW] Check for categorical_vars to preserve them and their dummies
      cat_vars <- attr(data, "categorical_vars")
      if (!is.null(cat_vars)) {
        for (cv in names(cat_vars)) {
          if (cv %in% available_vars) {
            # Add its dummy variables to available_vars so they aren't dropped
            available_vars <- unique(c(available_vars, cat_vars[[cv]]$dummies))
          }
        }
        # Keep only the subset that actually exists in data
        available_vars <- intersect(names(data), available_vars)
      }

      if (!quiet && length(data) > length(available_vars)) {
        message(sprintf(
          "Filtering data to %d relevant columns (out of %d) to optimize memory.",
          length(available_vars),
          length(data)
        ))
      }
      # Keep only relevant columns
      # If list, subset list. If df, subset df.
      if (is.data.frame(data)) {
        data <- data[, available_vars, drop = FALSE]
      } else {
        # [HIERARCHY FIX] If list data (hierarchical), do not subset the list by variable names.
        # Instead, filter columns within each dataframe in the list if needed.
        # For now, safe to keep the full list to avoid breaking level-to-level links.
      }

      # [NEW] Restore the attribute if it was present
      if (!is.null(cat_vars)) {
        attr(data, "categorical_vars") <- cat_vars
      }
    }
  }

  # Handle user-requested ordinal factors BEFORE categorical expansion
  if (!is.null(family)) {
    for (var in names(family)) {
      if (family[[var]] == "ordinal") {
        if (is.data.frame(data) && var %in% names(data)) {
          if (!is.ordered(data[[var]])) {
            data[[var]] <- factor(data[[var]], ordered = TRUE)
            if (!quiet) message(sprintf("Converted '%s' to ordered factor (family = 'ordinal')", var))
          }
        } else if (is.list(data) && !is.data.frame(data)) {
          for (i in seq_along(data)) {
            if (is.data.frame(data[[i]]) && var %in% names(data[[i]])) {
              if (!is.ordered(data[[i]][[var]])) {
                data[[i]][[var]] <- factor(data[[i]][[var]], ordered = TRUE)
                if (!quiet) message(sprintf("Converted '%s' to ordered factor (family = 'ordinal')", var))
              }
            }
          }
        }
      }
    }
  }

  # Save completely raw data before ANY categorical conversion happens
  # so that string/character tip labels remain intact for later PGLS matching
  original_raw_data_before_preprocess <- data

  # --- Automatic Data Cleaning (Handle Character/Factor Columns) ---
  data <- preprocess_categorical_vars(
    data,
    target_vars = model_vars,
    dummy_vars = fixed_predictors, # Categorical fixed predictors need dummies
    exclude_cols = id_col,
    quiet = quiet,
    expand_ordered = expand_ordered
  )

  # --- Hierarchical Data Detection & Validation ---
  # Data is hierarchical if it's a list (not dataframe)
  is_list_data <- is.list(data) && !is.data.frame(data)
  is_hierarchical <- FALSE
  hierarchical_info <- NULL

  if (is_list_data) {
    # Get all variables from fixed equations for auto-detection
    # We use fixed_eqs_temp to exclude random grouping variables from being 
    # assigned to a specific "home" level.
    eq_vars <- unique(unlist(lapply(fixed_eqs_temp, all.vars)))
    
    # Also explicitly exclude any variables used in random terms or global random arg
    all_random_groups <- unique(c(
      vapply(parsed_random_temp$random_terms, function(x) x$group, character(1)),
      global_random_vars
    ))
    eq_vars <- setdiff(eq_vars, all_random_groups)

    # Auto-detect if levels not provided
    if (is.null(levels)) {
      auto_result <- auto_detect_hierarchical(data, eq_vars, quiet = quiet)
      levels <- auto_result$levels

      # Use auto-detected multiscale hierarchy/link_vars if not explicitly provided
      if (is.null(multiscale)) {
        multiscale <- auto_result$hierarchy
      }
      if (is.null(link_vars)) {
        link_vars <- auto_result$link_vars
      }
    }

    # Now validate (with either provided or auto-detected values)
    if (!is.null(levels) && length(levels) > 0) {
      is_hierarchical <- TRUE
      
      # Handle alias
      if (is.null(multiscale) && !is.null(hierarchy)) multiscale <- hierarchy

      # [AUTO-DETECTION] Identify deterministic responses (LHS of I() equations)
      # Must happen BEFORE validation so the validator can skip them.
      det_responses <- character(0)
      for (eq in equations) {
          raw_eq <- formula(eq)
          if (length(raw_eq) >= 3) {
              resp <- as.character(raw_eq[[2]])
              rhs_str <- paste(deparse(raw_eq[[3]]), collapse = " ")
              if (grepl("^I\\(", trimws(rhs_str))) det_responses <- c(det_responses, resp)
          }
      }

      validate_hierarchical_data(
        data,
        levels,
        multiscale,
        link_vars,
        latent_vars = latent,
        equations = equations,
        deterministic_vars = det_responses
      )
      # Try to infer hierarchy from random effects if still not set
      if (is.null(multiscale)) {
        multiscale <- parse_hierarchy_from_random(random, data)
        if (is.null(multiscale)) {
          stop(
            "Multiscale data detected but 'multiscale' not specified. ",
            "Provide either:\n",
            "  1. 'multiscale' argument (e.g., \"site_year > individual\"), or\n",
            "  2. Nested random effects (e.g., ~(1|site/individual))"
          )
        }
      }


      # Store multiscale info for later use
      hierarchical_info <- list(
        data = data,
        levels = levels,
        hierarchy = multiscale,
        link_vars = link_vars,
        latent_vars = latent,
        deterministic_vars = det_responses # [NEW] Track deterministic responses
      )

      # Internal: Inject structural metadata if provided in sub-calls (e.g. from d-sep tests)
      if (!is.null(structure_multi)) {
        hierarchical_info$structure_multi <- structure_multi
      }
      if (!is.null(structure_levels)) {
        hierarchical_info$structure_levels <- structure_levels
      }

      if (!quiet) {
        message("Multiscale causal structure detected: ", multiscale)
      }
    }
  } else {
    # Data is not a data.frame (likely a list of dataframes or a standard JAGS list)
    if (!is.null(levels) && (!is.null(multiscale) || !is.null(hierarchy))) {
      is_hierarchical <- TRUE
      hierarchical_info <- list(
        data = data,
        levels = levels,
        hierarchy = if (!is.null(multiscale)) multiscale else hierarchy,
        link_vars = link_vars
      )
      
      # Internal: Inject structural metadata if provided in sub-calls
      if (!is.null(structure_multi)) {
        hierarchical_info$structure_multi <- structure_multi
      }
      if (!is.null(structure_levels)) {
        hierarchical_info$structure_levels <- structure_levels
      }
      
    }
  }

  # Capture row_ids early for hierarchical models, while 'data' is still the
  # raw input list of data frames (before preprocess_categorical_vars converts
  # character species columns to integer codes and before original_data is
  # overwritten at line ~1056).
  row_ids <- NULL
  if (is_hierarchical && !is.null(id_col) && is.list(data) && !is.data.frame(data)) {
    for (.lvl in names(data)) {
      .lvl_df <- data[[.lvl]]
      if (is.data.frame(.lvl_df) && id_col %in% names(.lvl_df)) {
        .col <- .lvl_df[[id_col]]
        if (is.character(.col) || is.factor(.col)) {
          row_ids <- as.character(.col)
          break
        }
      }
    }
    rm(.lvl, .lvl_df, .col)
  }

  if (is_hierarchical && !is.null(structure) && is.null(hierarchical_info$structure_levels)) {
      hierarchical_info$structure_levels <- auto_detect_structure_levels(structure, hierarchical_info, quiet = quiet)
  }

  # --- Random Effects Parsing ---
  # Extract (1|Group) and update equations to be fixed-effects only
  parsed_random <- extract_random_effects(equations)
  equations <- parsed_random$fixed_equations
  # Start with equation-specific random terms
  random_terms <- parsed_random$random_terms

  # Parse global random argument if provided
  if (!is.null(random)) {
    # --- DEPRECATION WARNING ---
    # Build the equivalent inline syntax to guide users on how to migrate.
    resp_vars <- vapply(
      parsed_random$fixed_equations,
      function(eq) as.character(eq[[2]]),
      character(1)
    )
    rand_char  <- as.character(random)
    rand_rhs   <- rand_char[length(rand_char)]
    rand_parts <- trimws(strsplit(rand_rhs, "\\+")[[1]])
    rand_parts <- rand_parts[grepl("\\|", rand_parts)]
    rand_inline <- paste(rand_parts, collapse = " + ")
    migration_lines <- vapply(
      resp_vars,
      function(r) paste0("  ", r, " ~ ... + ", rand_inline),
      character(1)
    )
    warning(
      "The 'random' argument is deprecated and will be removed in a future ",
      "version. Embed random effects directly in the equation formulas instead:\n",
      paste(migration_lines, collapse = "\n"), "\n",
      "See vignette('08_multiscale_models') for details.",
      call. = FALSE
    )

    global_random_terms <- parse_global_random(random, equations)
    # Combine with equation-specific terms
    random_terms <- c(random_terms, global_random_terms)

    # Deduplicate (avoid adding the same (1|Group) twice for the same response)
    if (length(random_terms) > 0) {
      keys <- vapply(
        random_terms,
        function(x) paste(x$response, x$group, sep = "|"),
        character(1)
      )
      random_terms <- random_terms[!duplicated(keys)]
    }
  }


  # --- Polynomial Term Extraction ---
  # Extract I(var^power) terms and expand formulas
  all_poly_terms <- get_all_polynomial_terms(equations)

  if (!is.null(all_poly_terms)) {
    # Expand formulas to replace I(x^2) with x_pow2
    equations <- lapply(equations, function(eq) {
      poly_terms <- extract_polynomial_terms(eq)
      expand_polynomial_formula(eq, poly_terms)
    })

    if (!quiet) {
      message(
        "Detected ",
        length(all_poly_terms),
        " polynomial term(s): ",
        paste(
          vapply(all_poly_terms, function(x) x$original, character(1)),
          collapse = ", "
        )
      )
    }

    # Auto-assign polynomial variables to the same level as their base variable
    if (is_hierarchical && !is.null(levels)) {
      for (poly in all_poly_terms) {
        base_var <- poly$base_var
        new_var <- poly$internal_name

        # Find level of base_var
        for (lvl_name in names(levels)) {
          if (base_var %in% levels[[lvl_name]]) {
            # Add the new poly var to this level
            levels[[lvl_name]] <- c(levels[[lvl_name]], new_var)
            break
          }
        }
      }

      # Update the info stored for later data retrieval
      hierarchical_info$levels <- levels
    }
  }
  # Initialize random structures (will be populated later)
  random_structures <- list()
  random_data_updates <- list()

  # Initialize result variables
  dsep_tests <- NULL
  dsep_results <- NULL
  dsep_correlations <- NULL

  # Handle global variability setting (e.g. variability = "reps")
  # If user provides a single string, apply it to all variables in equations
  if (
    !is.null(variability) &&
      is.character(variability) &&
      length(variability) == 1 &&
      is.null(names(variability))
  ) {
    global_type <- variability
    if (global_type %in% c("se", "reps")) {
      message(sprintf(
        "Global variability setting detected: applying '%s' to all variables.",
        global_type
      ))

      # Extract all variables from equations
      all_eq_vars <- unique(unlist(lapply(equations, all.vars)))

      # Exclude grouping variables from random effects
      grouping_vars <- character(0)
      if (!is.null(random)) {
        if (inherits(random, "formula")) {
          random_list <- list(random)
        } else {
          random_list <- random
        }

        for (r in random_list) {
          vars_in_random <- all.vars(r)
          grouping_vars <- c(grouping_vars, vars_in_random)
        }
      }

      # Also exclude Id col if provided
      if (!is.null(id_col)) {
        grouping_vars <- c(grouping_vars, id_col)
      }

      # Variables to apply variability to
      target_vars <- setdiff(all_eq_vars, grouping_vars)

      # Create named vector
      variability <- setNames(
        rep(global_type, length(target_vars)),
        target_vars
      )
    }
  }

  # --- Data Frame Preprocessing ---
  # If data is a data.frame, convert to list format expected by the model
  original_data <- data

  # --- Hierarchical Data Assembly ---
  # If hierarchical data provided, assemble full dataset for main model run
  if (is_hierarchical) {
    # Get all variables needed across all equations
    eq_vars <- unique(unlist(lapply(equations, all.vars)))

    # Add random effect grouping variables
    if (length(random_terms) > 0) {
      random_vars <- unique(vapply(
        random_terms,
        function(x) x$group,
        character(1)
      ))
      eq_vars <- unique(c(eq_vars, random_vars))
    }

    # Ensure all link_vars are in eq_vars. They define the hierarchy and may be
    # injected automatically as random effects during D-Separation tests.
    if (!is.null(hierarchical_info) && !is.null(hierarchical_info$link_vars)) {
      eq_vars <- unique(c(eq_vars, unlist(hierarchical_info$link_vars)))
    }

    # [NEW] Add categorical dummy variables to eq_vars so they aren't dropped
    # from the final hierarchical data list sent to JAGS.
    if (!is.null(attr(original_data, "categorical_vars"))) {
      cat_vars <- attr(original_data, "categorical_vars")
      for (cv in names(cat_vars)) {
        if (cv %in% eq_vars) {
          eq_vars <- unique(c(eq_vars, cat_vars[[cv]]$dummies))
        }
      }
    }

    # Remove latent variables (not in data)
    if (!is.null(latent)) {
      eq_vars <- setdiff(eq_vars, latent)
    }

    # Ensure base variables for polynomials are included (JAGS computes Age^2 from Age)
    # AND derived variables are excluded (so they aren't passed as data)
    if (!is.null(all_poly_terms)) {
      base_poly_vars <- vapply(
        all_poly_terms,
        function(x) x$base_var,
        character(1)
      )
      derived_poly_vars <- vapply(
        all_poly_terms,
        function(x) x$internal_name,
        character(1)
      )

      eq_vars <- unique(c(eq_vars, base_poly_vars))
      eq_vars <- setdiff(eq_vars, derived_poly_vars)
    }

    # Identify predictor variables (RHS) to generate dummies only for them
    rhs_vars <- unique(unlist(lapply(equations, function(eq) {
      if (length(eq) == 3) {
        all.vars(eq[[3]])
      } else {
        character(0)
      }
    })))

    # Save raw data with string/character values intact for d-sep PGLS tip matching
    original_raw_data <- original_raw_data_before_preprocess

    # Preprocess categorical variables in hierarchical data
    if (!is.null(hierarchical_info)) {
      hierarchical_info$data <- preprocess_categorical_vars(
        hierarchical_info$data,
        dummy_vars = rhs_vars,
        quiet = quiet
      )

      # Update 'data' variable if it's a list, so subsequent logic sees the attributes
      if (is.list(data) && !is.data.frame(data)) {
        data <- hierarchical_info$data
      }
    }

    # Include dummy variables in eq_vars so they are extracted by prepare_hierarchical_jags_data
    if (!is.null(attr(data, "categorical_vars"))) {
      cat_vars <- attr(data, "categorical_vars")
      # Extract RHS variables to identify which categorical vars are used as predictors
      rhs_vars <- unique(unlist(lapply(equations, function(eq) {
        if (length(eq) == 3) {
          all.vars(eq[[3]])
        } else {
          character(0)
        }
      })))

      for (cv_name in names(cat_vars)) {
        if (cv_name %in% rhs_vars) {
          # Add dummies to eq_vars
          eq_vars <- c(eq_vars, cat_vars[[cv_name]]$dummies)
        }
      }
      eq_vars <- unique(eq_vars)
    }

    # NEW PATH: Fully Hierarchical (Separate Loops)
    # Do NOT flatten. Prepare separate vectors and indices.

    prep_res <- prepare_hierarchical_jags_data(hierarchical_info, eq_vars)
    data <- prep_res$data_list

    # Restore categorical_vars attribute lost during list conversion
    if (exists("cat_vars") && !is.null(cat_vars)) {
      attr(data, "categorical_vars") <- cat_vars
    }

    # Add sample sizes to data list
    data <- c(data, prep_res$n_vec)

    if (!quiet) {
      message(
        paste0("Prepared hierarchical data for ", toupper(engine), ": "),
        paste(
          names(prep_res$n_vec),
          unlist(prep_res$n_vec),
          sep = "=",
          collapse = ", "
        )
      )
    }

    # original_data remains the raw un-indexed dataframes so tip labels match
    original_data <- original_raw_data
  }

  # --- Data Validation ---
  # Check that binomial response variables are 0/1 to avoid "Node inconsistent with parents" error
  if (!is.null(family)) {
    # If family is a single string but multiple equations, replicate it?
    # 'because' usually handles parsing family vector inside because_model, but we check here.
    # Assuming family corresponds to equations order or is named.

    # Simple check: scan equations and find their family
    # Note: 'family' argument handling in 'because' can be complex (vector vs single).
    # We'll use a simplified check using the names if possible, or position.

    responses <- vapply(
      equations,
      function(eq) as.character(formula(eq)[2]),
      character(1)
    )

    for (i in seq_along(responses)) {
      resp <- responses[i]
      fam <- NULL

      if (!is.null(names(family))) {
        if (resp %in% names(family)) {
          fam <- family[[resp]]
        }
      } else {
        if (length(family) == 1) {
          fam <- family
        } else if (length(family) == length(responses)) {
          # Ensure using double bracket if list, or single if vector
          if (is.list(family)) fam <- family[[i]] else fam <- family[i]
        }
      }

      if (!is.null(fam) && fam == "binomial" && resp %in% names(data)) {
        vals <- data[[resp]]
        if (any(!vals %in% c(0, 1, NA))) {
          unique_vals <- unique(vals[!is.na(vals)])
          stop(sprintf(
            "Binomial response variable '%s' contains invalid values: {%s}. Binomial variables must be strictly 0 or 1.",
            resp,
            paste(head(unique_vals, 5), collapse = ", ")
          ))
        }
      }
    }
  }

  # --- Random Effects Data Prep ---
  data <- prepare_random_effects_data(
    data = data, random_terms = random_terms, equations = equations,
    hierarchical_info = hierarchical_info, is_hierarchical = is_hierarchical,
    levels = levels, family = family, quiet = quiet
  )

  # --- Structure Processing ---
  struct_res <- process_because_structures(
    structure = structure, structures = structures, structure_obj = structure_obj,
    data = data, hierarchical_info = hierarchical_info,
    is_hierarchical = is_hierarchical, equations = equations,
    family = family, family_obj = family_obj, levels = levels,
    multiscale = multiscale, latent = latent, quiet = quiet
  )
  data       <- struct_res$data
  structures <- struct_res$structures

  # --- Structure Processing ---
  structures <- list()
  is_multiple <- FALSE
  N <- NULL

  # 1. Normalize Input to List
  if (is.null(structure)) {
    # Independent model - no structures to process
  } else if (is.matrix(structure)) {
    structures[["custom"]] <- structure
  } else if (is.list(structure) && !inherits(structure, "list")) {
    # It's an S3 object with list base - use class name
    class_name <- class(structure)[1]
    # Check if it's a multi-object type (contains multiple items)
    # Defense: Ignore list-based S3 objects that represent a single entity (like 'phylo' or 'spatial_knn')
    if (length(structure) > 1 && is.null(names(structure))) {
      is_multiple <- TRUE
      N_trees <- length(structure)
    }
    structures[[class_name]] <- structure
  } else if (
    is.list(structure) && (is.null(class(structure)) || identical(class(structure), "list"))
  ) {
    # Plain list of structures
    structures <- structure
    # Check for multi-objects in the list
    for (s in structures) {
      if (is.list(s) && length(s) > 1 && !is.matrix(s) && !inherits(s, "phylo") && !inherits(s, "because_structure")) {
        # Generic list with multiple items (likely replicates or multiPhylo)
        is_multiple <- TRUE
        if (!exists("N_trees")) N_trees <- length(structures)
      }
    }
  } else {
    # Any other S3 object (phylo, spatial_knn, etc.) - use class name
    class_name <- class(structure)[1]
    structures[[class_name]] <- structure
  }

  # Ensure N_trees is available if is_multiple is true
  if (is_multiple && !exists("N_trees")) {
     N_trees <- length(structures)
  }

  # Discover total N (number of observations) early
  if (is.null(N) && "N" %in% names(data)) {
    N <- if (is.list(data)) data$N[1] else data[["N"]][1]
  }
  if (is.null(N)) {
    # Try to get N from data - handle both data.frame and list cases
    if (is.data.frame(data) && nrow(data) > 0) {
      N <- nrow(data)
    } else if (is.list(data) && length(data) > 0) {
      # For hierarchical lists, N should be the row count of the FINEST grain level.
      first_obj <- data[[1]]
      if (is.data.frame(first_obj)) {
        N <- nrow(first_obj)
      } else if (is.vector(first_obj) || is.factor(first_obj)) {
        N <- length(first_obj)
      } else if (is.matrix(first_obj) || is.array(first_obj)) {
        N <- nrow(first_obj)
      }
    }
  }
  # [CLEANUP] Only include N if we are not in a complex hierarchical context where level-specific Ns take priority.
  # This silences the "Unused variable 'N' in data" warning in JAGS.
  if (is.null(hierarchical_info) && !"N" %in% names(data)) {
    data$N <- N
  }

  # [URGENT FIX] Pre-calculate level-specific counts so structure auto-detection works!
  # If we don't do this here, we can't match e.g. a 50x50 matrix to the 50-species level.
  if (!is.null(hierarchical_info)) {
    for (lvl_name in names(hierarchical_info$levels)) {
      z_name <- paste0("zeros_", lvl_name)
      if (is.null(data[[z_name]])) {
        data[[z_name]] <- rep(0, N)
      }
      n_name <- paste0("N_", lvl_name)
      if (is.null(data[[n_name]])) {
        if (is.data.frame(data) && !is.null(data[[lvl_name]])) {
          data[[n_name]] <- as.integer(length(unique(na.omit(data[[lvl_name]]))))[1]
        } else if (is.list(data) && !is.null(data[[lvl_name]])) {
          if (is.data.frame(data[[lvl_name]]) || is.matrix(data[[lvl_name]])) {
            data[[n_name]] <- as.integer(nrow(data[[lvl_name]]))[1]
          } else {
            data[[n_name]] <- as.integer(length(data[[lvl_name]]))[1]
          }
        } else {
          data[[n_name]] <- as.integer(N)[1]
        }
      }
      if (!is.null(data[[n_name]])) {
        data[[n_name]] <- as.integer(data[[n_name]])[1]
      }
      
      # Populate synonymous names (ID column names)
      lvl_vars <- hierarchical_info$levels[[lvl_name]]
      for (v_nm in lvl_vars) {
        v_title <- paste0(toupper(substring(v_nm, 1, 1)), substring(v_nm, 2))
        potential_names <- unique(c(v_nm, v_title, toupper(v_nm), paste0(v_nm, "ID"), paste0(v_title, "ID"), paste0(v_nm, "_id"), paste0(v_title, "_id")))
        for (p_nm in potential_names) {
            nn_name <- paste0("N_", p_nm)
            zz_name <- paste0("zeros_", p_nm)
            if (is.null(data[[nn_name]])) data[[nn_name]] <- data[[n_name]]
            if (is.null(data[[zz_name]])) data[[zz_name]] <- data[[z_name]]
        }
      }
    }
  }

  structure_names <- names(structures)
  if (is.null(structure_names) && length(structures) > 0) {
    structure_names <- paste0("Struct", seq_along(structures))
    names(structures) <- structure_names
  }

  # 2. Process Structures using S3 Generic
  if (length(structures) == 0) {
    # Independent Logic: Determine N from data
    if (is.null(N) || N == 0) {
      potential_objects <- Filter(
        function(x) is.vector(x) || is.factor(x) || is.matrix(x) || is.array(x),
        data
      )
      if (length(potential_objects) > 0) {
        obj <- potential_objects[[1]]
        N <- if (is.matrix(obj) || is.array(obj)) nrow(obj) else length(obj)
      }
    }
  } else {
    # --- HOTFIX: Inject raw string labels so prepare_structure_data can align tips ---
    # JAGS data is flattened to integers, which destroys tip labels. 
    # We temporarily inject them back into 'data' so they can be discovered.
    injected_raw_cols <- character(0)
    if (is_hierarchical && !is.null(hierarchical_info$structure_levels)) {
      for (s_name in names(hierarchical_info$structure_levels)) {
        lvl_name <- hierarchical_info$structure_levels[[s_name]]
        link_var <- if (!is.null(hierarchical_info$link_vars)) hierarchical_info$link_vars[[lvl_name]] else NULL
        if (!is.null(link_var) && is.data.frame(original_raw_data_before_preprocess[[lvl_name]])) {
           raw_col_name <- paste0(".raw_", link_var)
           data[[raw_col_name]] <- as.character(original_raw_data_before_preprocess[[lvl_name]][[link_var]])
           injected_raw_cols <- c(injected_raw_cols, raw_col_name)
        }
      }
    }

    structure_levels <- list()
    structure_multi <- list()
    # Use S3 Generic for Processing
    for (s_name in structure_names) {
      structure_obj <- structures[[s_name]]
      prep_res <- prepare_structure_data(structure_obj, data = data, optimize = TRUE, quiet = quiet, engine = engine, row_ids = row_ids)

      if (!is.null(prep_res$data_list)) {
        for (d_name in names(prep_res$data_list)) {
          custom_s_name <- get_structure_name_hook(structure_obj)
          prefixed_name <- if (d_name %in% c("Prec", "VCV", "multiVCV", custom_s_name)) paste0(d_name, "_", s_name) else d_name
          
          # --- HOTFIX: Prevent Python/JAGS crashes ---
          # prepare_structure_data might return character/factor vectors (e.g. aligned tip labels).
          # We MUST NOT let these overwrite the integer index arrays (e.g. data$Species)
          # nor be added to the data list, as NumPyro strictly requires numeric arrays.
          if (is.character(prep_res$data_list[[d_name]]) || is.factor(prep_res$data_list[[d_name]])) {
             next
          }
          
          data[[prefixed_name]] <- prep_res$data_list[[d_name]]
        }
      } else {
        structures[[s_name]] <- NULL
        structure_names <- setdiff(structure_names, s_name)
        next
      }

      current_N <- NULL
      is_this_one_multi <- FALSE
      for (obj in prep_res$data_list) {
        if (is.matrix(obj) && nrow(obj) == ncol(obj)) {
          current_N <- nrow(obj)
          break
        } else if (is.array(obj) && length(dim(obj)) == 3) {
          dims <- dim(obj)
          if (dims[1] == dims[2]) {
             current_N <- dims[1]
             is_this_one_multi <- TRUE
             if (is_multiple && !exists("N_trees")) N_trees <- dims[3]
             break
          } else if (dims[2] == dims[3]) {
             current_N <- dims[2]
             is_this_one_multi <- TRUE
             if (is_multiple && !exists("N_trees")) N_trees <- dims[1]
             break
          }
        }
      }

      if (is.null(current_N)) {
        n_attr <- attr(structure_obj, "n")
        if (is.numeric(n_attr)) current_N <- n_attr
      }

      # [FIX] Robust Level Matching
      # Prioritize name-based matching (e.g. if the structure in the tree list is named "Species")
      s_level <- NULL
      if (!is.null(hierarchical_info) && !is.null(current_N)) {
        # 1. Try Name Matching first
        if (s_name %in% names(hierarchical_info$levels)) {
            s_level <- s_name
        } else {
            # 2. Fallback to Dimension Matching
            for (lvl in names(hierarchical_info$levels)) {
              n_name <- paste0("N_", lvl)
              if (!is.null(data[[n_name]]) && data[[n_name]] == current_N) {
                s_level <- lvl
                break
              }
            }
        }
      }
      structure_levels[[s_name]] <- s_level
      structure_multi[[s_name]]  <- is_this_one_multi

      if (!is.null(current_N)) {
        if (is.null(N) || N == 0) {
          N <- current_N
        } else if (N != current_N && is.null(hierarchical_info)) {
          stop(paste("Dimension mismatch in structure:", s_name))
        }
      }

      if (!is.null(prep_res$structure_object)) {
        structures[[s_name]] <- prep_res$structure_object
      }
    }

    if (!is.null(hierarchical_info)) {
      # Merge: prefer the already name-matched auto-detect result for each s_name,
      # only overwriting with the newly computed (dimension-based) value when the
      # auto-detect didn't produce a result for that structure.
      existing_sl <- hierarchical_info$structure_levels
      for (s_name in names(structure_levels)) {
        # Only override if auto-detect has no entry or if dimension match is more specific
        if (is.null(existing_sl[[s_name]])) {
          existing_sl[[s_name]] <- structure_levels[[s_name]]
        }
        # If both give a result, keep the existing (name-based) auto-detect result
      }
      hierarchical_info$structure_levels <- existing_sl
      hierarchical_info$structure_multi  <- structure_multi
    }
  }

  # --- [NEW] Map Link Variables (IDs) to Counts ---
  # Ensures JAGS finds loop bounds for grouping variables not listed in 'levels'
  if (!is.null(hierarchical_info$link_vars) && !is.null(hierarchical_info$data)) {
    for (lk_var in hierarchical_info$link_vars) {
      # Find which level (dataframe) contains this link variable
      # If multiple, find the one with the fewest rows (the level where it's defined)
      potential_lvls <- character(0)
      for (l_nm in names(hierarchical_info$data)) {
         if (lk_var %in% colnames(hierarchical_info$data[[l_nm]])) {
            potential_lvls <- c(potential_lvls, l_nm)
         }
      }
      
      if (length(potential_lvls) > 0) {
        # Pick the level with minimum rows
        r_counts <- vapply(potential_lvls, function(l) nrow(hierarchical_info$data[[l]]), numeric(1))
        best_lvl <- potential_lvls[which.min(r_counts)]
        best_n <- r_counts[best_lvl]
        
        # Generate names (Raw, Title, Upper)
        v_title <- paste0(toupper(substring(lk_var, 1, 1)), substring(lk_var, 2))
        pot_names <- unique(c(lk_var, v_title, toupper(lk_var)))
        
        for (p_nm in pot_names) {
            nn_name <- paste0("N_", p_nm)
            if (is.null(data[[nn_name]])) {
              data[[nn_name]] <- best_n
            }
        }
      }
    }
  }

  # Only add 'zeros' vector if using ZIP or ZINB (Poisson trick)
  # Only add 'zeros' vector if using ZIP or ZINB (Poisson trick)
  # OR if we have structures (which often use dmnorm(zeros, ...))
  if (!is.null(N)) {
    # Check if 'zeros' already exists (use exact name match)
    has_zeros <- "zeros" %in% names(data)
    if (!has_zeros) {
      needs_zeros <- length(structures) > 0 ||
        any(vapply(
          names(family),
          function(v) {
            # Dispatch on the specific family object for this variable
            fam_name <- family[[v]]
            fam_obj_v <- get_family(fam_name)
            needs_zero_inflation_hook(fam_obj_v, v)
          },
          logical(1)
        ))

      if (needs_zeros) {
        data[["zeros"]] <- rep(0, N)
      }
    }
  }

  # ID is currently unused in model templates, removing to avoid JAGS warnings
  # ID2 is used for pairwise induced correlations (Wishart priors)
  data$ID2 <- diag(2)

  # Handle multinomial and ordinal data
  if (!is.null(family)) {
    # If family is provided but unnamed, try to auto-assign if there is only one response
    if (is.null(names(family))) {
      response_vars <- unique(vapply(
        equations,
        function(eq) as.character(all.vars(eq[[2]])[1]),
        character(1)
      ))
      if (length(family) == 1 && length(response_vars) == 1) {
        names(family) <- response_vars
        message(sprintf(
          "Auto-assigned family '%s' to response variable '%s'",
          family,
          response_vars
        ))
      } else if (length(family) == length(response_vars)) {
        # Riskier, but if lengths match, assume order
        # Better to warn and ask for names.
        warning(
          "Argument 'family' is unnamed. Please provide a named vector like c(Response = 'binomial'). Assuming defaults (Gaussian) for safety."
        )
      } else {
        warning(
          "Argument 'family' is unnamed and length does not match response variables. Ignoring."
        )
      }
    }

    # Auto-fix residual variance for non-Gaussian distributions if not specified
    for (var in names(family)) {
      dist <- family[[var]]
      if (dist %in% c("binomial", "multinomial", "ordinal")) {
        should_fix <- FALSE
        if (is.null(fix_residual_variance)) {
          should_fix <- TRUE
          fix_residual_variance <- c()
        } else if (
          is.numeric(fix_residual_variance) &&
            length(fix_residual_variance) == 1 &&
            is.null(names(fix_residual_variance))
        ) {
          # It's a global fix (e.g. fix=1), so it applies to this var too. No action needed.
          should_fix <- FALSE
        } else if (!var %in% names(fix_residual_variance)) {
          should_fix <- TRUE
        }

        if (should_fix) {
          # Append to fixed variance vector
          new_fix <- setNames(1, var)
          fix_residual_variance <- c(fix_residual_variance, new_fix)

          if (!quiet) {
            message(sprintf(
              "Note: Fixing residual variance of '%s' (%s) to 1 for identifiability.",
              var,
              dist
            ))
          }
        }
      }
    }

    # Identify categorical variables and their K levels
    cat_vars_metadata <- attr(data, "categorical_vars")

    for (var in names(family)) {
      if (family[[var]] %in% c("multinomial", "ordinal")) {
        if (!var %in% names(data)) {
          stop(paste(
            family[[var]],
            "variable",
            var,
            "not found in data."
          ))
        }

        # Determine K (number of levels)
        # PRIORITIZE existing categorical metadata if available
        if (!is.null(cat_vars_metadata) && var %in% names(cat_vars_metadata)) {
          K <- length(cat_vars_metadata[[var]]$levels)
          val <- data[[var]]
          if (is.factor(val)) {
             data[[var]] <- as.integer(val)
          } else if (is.character(val)) {
             data[[var]] <- as.integer(factor(val, levels = cat_vars_metadata[[var]]$levels))
          }
        } else {
          # Fallback: Detect from current data vector
          val <- data[[var]]
          if (is.factor(val)) {
            K <- nlevels(val)
            data[[var]] <- as.integer(val)
          } else {
            # Assume it's already integer or character
            val <- as.factor(val)
            K <- nlevels(val)
            data[[var]] <- as.integer(val)
          }
        }

        if (family[[var]] == "multinomial" && K < 3) {
          warning(paste(
            "Multinomial variable",
            var,
            "has fewer than 3 levels. Consider using binomial."
          ))
        }

        if (family[[var]] == "ordinal" && K < 3) {
          warning(paste(
            "Ordinal variable",
            var,
            "has fewer than 3 levels. Consider using binomial."
          ))
        }

        # Pass K to JAGS (auto-detected from data)
        K_name <- paste0("K_", var)
        if (!K_name %in% names(data)) {
          data[[K_name]] <- K
          if (!quiet) {
            message(sprintf(
              "Auto-detected K_%s = %d from %s variable '%s'",
              var,
              K,
              family[[var]],
              var
            ))
          }
        }
      }
    }
  }

  if (!is.null(family)) {
    for (var_name in names(family)) {
      dist_type <- family[[var_name]]
      if (dist_type %in% c("zip", "zinb")) {
        zeros_name <- paste0("zeros_", var_name)
        if (is.null(data[[zeros_name]])) {
          data[[zeros_name]] <- rep(0, N)
        }
      }
    }
  }

  # Check for missing data
  all_vars <- unique(unlist(lapply(equations, function(eq) {
    c(all.vars(eq[[3]]), all.vars(eq[[2]]))
  })))

  # [FIX] Ensure parent categorical variables are included in all_vars
  # so their missing values (NAs) are detected and handled even if they
  # were expanded into dummy variables in the formulas.
  cat_metadata <- attr(data, "categorical_vars")
  if (!is.null(cat_metadata)) {
    for (parent in names(cat_metadata)) {
      dummies <- cat_metadata[[parent]]$dummies
      if (any(dummies %in% all_vars)) {
        all_vars <- unique(c(all_vars, parent))
      }
    }
  }

  response_vars <- unique(vapply(
    equations,
    function(eq) as.character(all.vars(eq[[2]])[1]),
    character(1)
  ))
  predictor_only_vars <- setdiff(all_vars, response_vars)

  # Detect variables with missing data
  response_vars_with_na <- character(0)
  predictor_vars_with_na <- character(0)

  for (var in all_vars) {
    # [MULTISCALE FIX] Check inside list of datasets if data is multiscale
    var_data <- if (var %in% names(data)) {
      data[[var]]
    } else if (is.list(data) && !is.data.frame(data)) {
      # Search through levels to find the variable
      found_val <- NULL
      for (lvl_name in names(data)) {
        if (is.data.frame(data[[lvl_name]]) && var %in% colnames(data[[lvl_name]])) {
          found_val <- data[[lvl_name]][[var]]
          break
        }
      }
      found_val
    } else {
      NULL
    }

    if (!is.null(var_data)) {
      if (
        !is.matrix(var_data) && any(is.na(var_data)) && !all(is.na(var_data))
      ) {
        if (var %in% response_vars) {
          response_vars_with_na <- c(response_vars_with_na, var)
        } else {
          predictor_vars_with_na <- c(predictor_vars_with_na, var)
        }
      }
    }
  }

  # Handle predictor-only variables with missing data
  if (length(predictor_vars_with_na) > 0) {
    if (!quiet) {
      message(
        "Note: Detected missing data in predictor-only variables: ",
        paste(predictor_vars_with_na, collapse = ", "),
        "\nAutomatically adding intercept-only equations (e.g., ",
        predictor_vars_with_na[1],
        " ~ 1) to enable imputation."
      )
    }

    # Add intercept-only equations
    for (var in predictor_vars_with_na) {
      new_eq <- stats::as.formula(paste(var, "~ 1"))
      equations <- c(equations, list(new_eq))
    }

    # Treat them as responses now
    response_vars_with_na <- c(response_vars_with_na, predictor_vars_with_na)
  }

  # For variables with missing data, we'll use the GLMM (Latent Variable) approach

  # Auto-detect variability from data column names (user-friendly)
  # Look for patterns: X_se, X_obs or matrix columns
  auto_variability <- list()

  for (var in all_vars) {
    # Skip if already in manual variability specification
    if (!is.null(variability) && var %in% c(names(variability), variability)) {
      next
    }

    # Extension Hook: Variability priority handling
    if (!is.null(family_obj)) {
      v_type <- get_variability_type_hook(family_obj, var)
      if (!is.null(v_type)) {
        auto_variability[[var]] <- v_type
        if (!quiet) {
          message(sprintf(
            "Extension-detected: '%s' is a specialized family -> using '%s' mode.",
            var,
            v_type
          ))
        }
        next
      }
    }

    # Check for SE pattern (X_se)
    se_name <- paste0(var, "_se")
    sd_name <- paste0(var, "_sd")

    if (se_name %in% names(data)) {
      auto_variability[[var]] <- "se"
      if (!quiet) {
        message(sprintf(
          "Auto-detected: '%s' has standard errors in '%s'",
          var,
          se_name
        ))
      }

      # Check for repeated measures pattern (X_obs or matrix)
      obs_name <- paste0(var, "_obs")
      if (var %in% names(data) && is.matrix(data[[var]])) {
        auto_variability[[var]] <- "reps"
        if (!quiet) {
          message(sprintf(
            "Auto-detected: '%s' has repeated measures (matrix format)",
            var
          ))
        }
      } else if (obs_name %in% names(data)) {
        auto_variability[[var]] <- "reps"
        if (!quiet) {
          message(sprintf(
            "Auto-detected: '%s' has repeated measures in '%s'",
            var,
            obs_name
          ))
        }
      }
    }
  }

  # Merge auto-detected with manual specification (manual takes precedence)
  if (length(auto_variability) > 0) {
    if (is.null(variability)) {
      variability <- auto_variability
    } else {
      # Convert variability to named list if needed
      if (is.null(names(variability))) {
        # This case (unnamed vector but length > 1) is ambiguous or unsupported by global logic
        # We'll treat as "names missing" warning or error?
        # For backward compatibility / safety, just set names to values?
        # Actually, existing code: variability <- setNames(rep(NA, ...)) seems wrong if we passed values like c("reps", "se")
        # But previous logic (line 958 in original) did: variability <- setNames(rep(NA, length(variability)), variability)
        # This implied variability was a vector of NAMES.
        # But wait, the doc says variability is "c(Var = 'type')".
        # If unnamed, `variability` contains TYPES? Or NAMES?
        # Original Doc line 11 (approx): "If unnamed, it defaults to 'se' for all specified variables."
        # Meaning: variability = c("Var1", "Var2") -> Var1="se", Var2="se".
        # My new global logic handles length==1 separately.
        # If length > 1 and unnamed, we assume standard behavior (list of vars, default type "se")
        variability <- setNames(rep("se", length(variability)), variability)
      }

      # Merge: manual overrides auto
      for (var in names(auto_variability)) {
        if (!var %in% names(variability)) {
          variability[[var]] <- auto_variability[[var]]
        }
      }
    }
  }
  # Handle variability data
  variability_list <- list()
  if (!is.null(variability)) {
    for (var_name in names(variability)) {
      var_spec <- variability[[var_name]]

      # Parse specification: can be "se"/"reps" or list(type="se", se_col="X_SD")
      if (is.list(var_spec)) {
        # Extended format with custom column names
        type <- var_spec$type
        custom_se_col <- var_spec$se_col
        custom_obs_col <- var_spec$obs_col
        custom_mean_col <- var_spec$mean_col
      } else {
        # Simple format: just the type
        type <- as.character(var_spec)
        custom_se_col <- NULL
        custom_obs_col <- NULL
        custom_mean_col <- NULL
      }

      # Validate type
      if (!type %in% c("se", "reps")) {
        stop(paste(
          "Invalid variability type for",
          var_name,
          "- must be 'se' or 'reps', got:",
          type
        ))
      }

      variability_list[[var_name]] <- type

      if (type == "se") {
        # Determine SE column name (custom or standard)
        se_col <- custom_se_col %||% paste0(var_name, "_se")
        mean_col <- custom_mean_col %||% paste0(var_name, "_mean")

        # Check if SE column exists
        if (!se_col %in% names(data)) {
          stop(paste(
            "Variable",
            var_name,
            "specified as 'se' type but column",
            se_col,
            "not found in data."
          ))
        }

        # Handle mean column
        if (!mean_col %in% names(data)) {
          if (var_name %in% names(data)) {
            # Rename var to var_mean
            data[[mean_col]] <- data[[var_name]]
            data[[var_name]] <- NULL
          } else {
            stop(paste(
              "Variable",
              var_name,
              "specified as 'se' type but neither",
              var_name,
              "nor",
              mean_col,
              "found in data."
            ))
          }
        } else {
          # If both exist, ensure var is removed (it's a latent parameter now)
          if (var_name %in% names(data)) data[[var_name]] <- NULL
        }

        # Rename custom column to standard name if needed
        if (se_col != paste0(var_name, "_se")) {
          data[[paste0(var_name, "_se")]] <- data[[se_col]]
        }
      } else if (type == "reps") {
        # Determine obs column name (custom or standard)
        obs_col <- custom_obs_col %||% paste0(var_name, "_obs")
        nrep_name <- paste0("N_reps_", var_name)

        # Check if obs column exists
        if (!obs_col %in% names(data)) {
          # Also allow data.frame (as user detection histories often come as DF)
          # Case 1: Matrix or Data Frame (Single Species or explicit multi-dimensional)
          if (
            var_name %in%
              names(data) &&
              (is.matrix(data[[var_name]]) || is.data.frame(data[[var_name]]))
          ) {
            data[[obs_col]] <- as.matrix(data[[var_name]])
            data[[var_name]] <- NULL
          } else if (
            var_name %in%
              names(data) &&
              is.list(data[[var_name]]) &&
              !is.data.frame(data[[var_name]])
          ) {
            # Case 2: List of Matrices (Multispecies Bundle)
            # Verify elements are matrices/dfs
            elem_valid <- all(sapply(data[[var_name]], function(x) {
              is.matrix(x) || is.data.frame(x)
            }))
            if (!elem_valid) {
              stop(
                "Elements of '",
                var_name,
                "' list must be matrices or data frames."
              )
            }

            # Convert to 3D Array [Sites, Reps, Species]
            # Assumption: All matrices have same dimensions
            tryCatch(
              {
                # Convert all to matrices first
                mat_list <- lapply(data[[var_name]], as.matrix)
                arr_3d <- simplify2array(mat_list)
                data[[obs_col]] <- arr_3d
                data[[var_name]] <- NULL

                if (!quiet) {
                  message(
                    "Converted list of matrices '",
                    var_name,
                    "' to 3D array (",
                    paste(dim(arr_3d), collapse = "x"),
                    ")."
                  )
                }
              },
              error = function(e) {
                stop(
                  "Failed to convert list '",
                  var_name,
                  "' to 3D array. Ensure all matrices have identical dimensions. Error: ",
                  e$message
                )
              }
            )
          } else {
            stop(paste(
              "Variable",
              var_name,
              "specified as 'reps' type but column",
              obs_col,
              "(as matrix) not found in data."
            ))
          }
        }

        # Rename custom column to standard name if needed
        if (obs_col != paste0(var_name, "_obs")) {
          data[[paste0(var_name, "_obs")]] <- data[[obs_col]]
        }

        # Calculate N_reps if not provided
        if (!nrep_name %in% names(data)) {
          mat <- data[[paste0(var_name, "_obs")]]
          # Count non-NA values per row
          n_reps <- apply(mat, 1, function(x) sum(!is.na(x)))
          data[[nrep_name]] <- n_reps

          # Compact matrix (move non-NAs to left) to ensure 1:N_reps indexing works
          compact_mat <- matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
          for (i in seq_len(nrow(mat))) {
            vals <- mat[i, !is.na(mat[i, ])]
            if (length(vals) > 0) {
              compact_mat[i, seq_along(vals)] <- vals
            }
          }
          data[[paste0(var_name, "_obs")]] <- compact_mat
        }
      }

      # Ensure var is removed (latent)
      if (var_name %in% names(data)) data[[var_name]] <- NULL
    }
  }

    # Auto-detect latent variables: variables in equations but not in data
    if (is.null(latent)) {
      vars_in_equations <- unique(unlist(lapply(equations, all.vars)))
      vars_in_data <- names(data)

      # Find variables that appear in equations but not in data
      potential_latents <- setdiff(vars_in_equations, vars_in_data)

      # Extension Hook: Remove specialized variables from potential latents
      potential_latents <- dsep_potential_latent_hook(
        family_obj,
        potential_latents
      )

      # Exclude polynomial internal variables (they're deterministic, not latent)
      if (!is.null(all_poly_terms)) {
        poly_internal_names <- sapply(all_poly_terms, function(x) {
          x$internal_name
        })
        potential_latents <- setdiff(potential_latents, poly_internal_names)
      }

      # Exclude categorical parent variables (they are replaced by dummies but are still in the model equations)
      cat_metadata <- attr(data, "categorical_vars")
      if (!is.null(cat_metadata)) {
        potential_latents <- setdiff(potential_latents, names(cat_metadata))
      }

      if (length(potential_latents) > 0) {
        # Auto-detect latent variables
        latent <- potential_latents

        if (!quiet) {
          msg <- paste0(
            "Auto-detected latent variable(s): ",
            paste(latent, collapse = ", "),
            "\n(Variables in equations but not in data will be treated as latent.)"
          )
          if (dsep) {
            msg <- paste0(msg, "\nGenerating m-separation tests for MAG...")
          }
          message(msg)
        }
      }
    }

    if (!quiet && dsep) {
      if (!is.null(latent)) {
        message(
          "Generating m-separation tests (MAG with latent variables)..."
        )
      } else {
        message("Generating d-separation tests...")
      }
    }

  # --- D-Separation Tests ---
  dsep_res <- run_because_dsep(
    dsep = dsep, engine = engine, equations = equations,
    data = data, family = family, structures = structures,
    structure_obj = structure_obj,
    hierarchical_info = hierarchical_info, is_hierarchical = is_hierarchical,
    random_terms = random_terms, levels = levels, multiscale = multiscale,
    link_vars = link_vars, latent = latent, latent_method = latent_method,
    id_col = id_col, variability = variability,
    all_poly_terms = all_poly_terms,
    fixed_equations_temp = fixed_equations_temp,
    induced_cors_in = induced_cors, dsep_max_obs = dsep_max_obs,
    aggregate_crossscale = aggregate_crossscale,
    parallel = parallel, n.cores = n.cores, n.chains = n.chains,
    n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.adapt = n.adapt,
    ic_recompile = ic_recompile, fix_residual_variance = fix_residual_variance,
    quiet = quiet, priors = priors, monitor_mode = monitor_mode,
    expand_ordered = expand_ordered, nimble_samplers = nimble_samplers,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    prior_scale_fixed = prior_scale_fixed, verbose = verbose
  )
  dsep_tests <- dsep_res$dsep_tests
  induced_cors <- dsep_res$induced_cors
  equations    <- dsep_res$equations

  # --- Ensure zero_vec and ID2 are in data if needed (for multivariate priors) ---
  # Check if we need zero_vec (matches logic in because_model.R)
  need_zero_vec_data <- (!is.null(induced_cors) && length(induced_cors) > 0) || 
                        (!is.null(structures) && length(structures) > 0) ||
                        (!is.null(family) && (any(family == "multinomial") || any(family == "ordinal")))

  if (need_zero_vec_data) {
    # Check if already provided (e.g. by hierarchical_prep), otherwise add it
    if (is.null(data$zero_vec)) {
      # Determine a safe max length
      # Use 1000 as a safe minimum, or actual N if available
      n_max <- 1000
      if (!is.null(data$N)) n_max <- max(n_max, data$N)
      # In hierarchical models, check for level sample sizes
      n_names <- grep("^N_", names(data), value = TRUE)
      if (length(n_names) > 0) {
        for (nn in n_names) n_max <- max(n_max, data[[nn]])
      }
      
      data$zero_vec <- rep(0, n_max)
    }
    
    # ID2 is used for Wishart priors in induced correlations
    if (is.null(data$ID2)) {
      data$ID2 <- diag(2)
    }
  }

  # --- Normalize Family Map for Metadata ---
  # Ensure family is a named character vector containing ALL response variables.
  # This prevents "subscript out of bounds" in downstream S3 methods (pp_check, predict).
  temp_responses <- vapply(equations, function(eq) as.character(formula(eq)[2]), character(1))
  
  if (is.null(family)) {
    family <- setNames(rep("gaussian", length(temp_responses)), temp_responses)
  } else {
    # If family is a single string or was partially named, expand it
    if (is.null(names(family)) && length(family) == 1) {
       family <- setNames(rep(family, length(temp_responses)), temp_responses)
    } else {
       # Ensure every variable is present
       for (r in temp_responses) {
         if (!(r %in% names(family))) {
           family[[r]] <- "gaussian"
         }
       }
    }
  }

  # Deduplicate random structures: prevent names in 'structure' from being 
  # treated as generic random intercepts by 'because_model'
  r_names <- names(random_structures)
  if (!is.null(structures)) {
    # Case-insensitive setdiff
    struct_names_lower <- tolower(names(structures))
    keep_idx <- !tolower(r_names) %in% struct_names_lower
    r_names <- r_names[keep_idx]
  }

  # JAGS model code
  model_output <- because_model(
    equations = equations,
    is_multi_structure = is_multiple,
    variability = variability_list,
    family = if (!is.null(family)) as.list(family) else NULL,
    vars_with_na = response_vars_with_na,
    induced_correlations = induced_cors,
    latent = latent,
    standardize_latent = standardize_latent,
    fix_latent = fix_latent,
    structures = structures,
    random_structure_names = r_names,
    random_terms = random_terms,
    poly_terms = all_poly_terms,
    categorical_vars = if (!is.null(attr(data, "categorical_vars"))) {
      attr(data, "categorical_vars")
    } else {
      NULL
    },
    priors = priors,
    hierarchical_info = if (is_hierarchical) hierarchical_info else NULL,
    engine = engine,
    quiet = quiet
  )

  model_string <- model_output$model
  parameter_map <- model_output$parameter_map
  if (engine == "numpyro") {
    eq_strings <- sapply(equations, function(eq) paste(deparse(eq), collapse=" "))
    
    # Process structures for NumPyro
    py_structures <- list()
    for (s_name in names(structures)) {
      s_obj <- structures[[s_name]]
      
      # Determine matrix name in the prepared data
      custom_s_name <- get_structure_name_hook(s_obj)
      mat_name <- if (!is.null(custom_s_name)) paste0(custom_s_name, "_", s_name) else paste0("VCV_", s_name)
      if (is.null(data[[mat_name]]) && !is.null(data[[paste0("Prec_", s_name)]])) {
        mat_name <- paste0("Prec_", s_name)
      }
      
      if (!is.null(data[[mat_name]])) {
        # Ask the extension package for the Python JAX code
        py_code <- numpyro_structure_definition(s_obj, engine = "numpyro")
        if (!is.null(py_code)) {
          # Compile into a Python function using reticulate
          env <- reticulate::py_run_string(py_code)
          funcs <- names(env)
          expected_name <- paste0(s_name, "_transform")
          func_name <- if (expected_name %in% funcs) expected_name else funcs[!funcs %in% c("numpyro", "jnp", "jax", "dist", "np", "r")][1]
          if (!is.null(func_name) && (func_name %in% names(env))) {
            py_structures[[s_name]] <- list(
              matrix = data[[mat_name]],
              transform_func = env[[func_name]],
              type = class(s_obj)[1]
            )
          }
        } else {
          # Fallback to pure matrix
          py_structures[[s_name]] <- data[[mat_name]]
        }
        
        # Append + (1 | s_name) to valid endogenous equations
        for (i in seq_along(equations)) {
           eq <- equations[[i]]
           response <- trimws(strsplit(deparse(eq), "~")[[1]][1])
           
           # Apply if dimension matches or valid level
           is_valid <- TRUE
           if (exists("is_hierarchical") && is_hierarchical && exists("hierarchical_info")) {
               s_lvl <- hierarchical_info$structure_levels[[s_name]]
               tryCatch({
                   resp_lvl <- infer_variable_level(response, hierarchical_info$levels, data = NULL, equations = equations, latent = latent, hierarchy = hierarchical_info$hierarchy)
                   if (is.null(resp_lvl)) {
                       is_valid <- FALSE
                   } else if (!is_valid_structure_mapping_dsep(s_lvl, resp_lvl, hierarchical_info)) {
                       is_valid <- FALSE
                   }
               }, error = function(e) {
                   is_valid <<- FALSE   # <<- assigns to enclosing loop scope, not local handler scope
               })
           }
           
           if (is_valid) {
               eq_str <- eq_strings[i]
               if (!grepl(paste0("\\(1\\s*\\|\\s*", s_name, "\\)"), eq_str)) {
                   eq_strings[i] <- paste0(eq_str, " + (1|", s_name, ")")
               }
           }
        }
      }
    }
    
    flat_data <- flatten_for_python(data)
    
    # Ensure zero-indexing for Python categorical variables
    for (s_name in names(structures)) {
        # If indexing array is missing (e.g. 1-to-1 row mapping like JAGS), generate it!
        if (is.null(flat_data[[s_name]])) {
            N_val <- NULL
            if (exists("is_hierarchical") && is_hierarchical && exists("hierarchical_info")) {
                s_lvl <- hierarchical_info$structure_levels[[s_name]]
                if (!is.null(s_lvl)) {
                    N_val <- flat_data[[paste0("N_", s_lvl)]]
                }
            }
            if (is.null(N_val)) {
                N_val <- if (!is.null(flat_data[["N"]])) flat_data[["N"]] else length(flat_data[[1]])
            }
            flat_data[[s_name]] <- 0:(N_val - 1)
        }
    }
    idx_vars <- grep("_idx", names(flat_data), value = TRUE)
    idx_vars <- c(idx_vars, names(structures))
    for (eq in equations) {
        eq_str <- if (is.character(eq)) eq else paste(deparse(eq), collapse=" ")
        matches <- regmatches(eq_str, gregexpr("\\(1\\s*\\|\\s*[^)]+\\)", eq_str))[[1]]
        for (m in matches) {
            grp <- trimws(strsplit(m, "\\|")[[1]][2])
            grp <- gsub("\\)", "", grp)
            idx_vars <- c(idx_vars, grp)
        }
    }
    idx_vars <- unique(idx_vars)
    for (s_name in idx_vars) {
        if (s_name %in% names(flat_data) && min(flat_data[[s_name]], na.rm=TRUE) >= 1) {
            flat_data[[s_name]] <- as.integer(flat_data[[s_name]] - 1L)
        }
    }
    
    # Re-inject standard random terms that were stripped for JAGS
    if (length(random_terms) > 0) {
        for (rt in random_terms) {
            for (i in seq_along(equations)) {
                resp <- trimws(strsplit(deparse(equations[[i]]), "~")[[1]][1])
                if (resp == rt$response) {
                    re_str <- paste0("\\(1\\s*\\|\\s*", rt$group, "\\)")
                    if (!grepl(re_str, eq_strings[i])) {
                        eq_strings[i] <- paste0(eq_strings[i], " + (1|", rt$group, ")")
                    }
                }
            }
        }
    }

    py_result <- run_numpyro_model(
      eq_strings = eq_strings,
      flat_data = flat_data,
      family = if (!is.null(family)) as.list(family) else NULL,
      priors = NULL,
      py_structures = py_structures,
      n_chains = n.chains,
      n_iter = n.iter - n.burnin,
      n_warmup = n.burnin,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      prior_scale_fixed = prior_scale_fixed,
      quiet = quiet
    )
    mcmc_samples <- format_numpyro_samples(py_result)
    
    result <- list(
      equations = equations,
      model      = NULL,        # No live model object for NumPyro
      model_code = model_string, # JAGS-equivalent string (kept for reference)
      numpyro_code = if (!is.null(py_result$model_code)) py_result$model_code else NULL,
      parameter_map = parameter_map,
      samples = mcmc_samples,
      data = data,
      dsep = NULL,
      priors = priors,
      hierarchical_info = if (is_hierarchical) hierarchical_info else NULL,
      engine = engine,
      quiet = quiet
    )
    if (WAIC && !is.null(py_result$waic)) {
      waic_df <- data.frame(
        Estimate = c(py_result$waic$elpd_waic$Estimate, py_result$waic$p_waic$Estimate, py_result$waic$waic$Estimate),
        SE = c(py_result$waic$elpd_waic$SE, py_result$waic$p_waic$SE, py_result$waic$waic$SE),
        row.names = c("elpd_waic", "p_waic", "waic")
      )
      
      attr(waic_df, "pointwise") <- data.frame(
        elpd_waic = as.numeric(py_result$waic$pointwise$elpd_waic_i),
        p_waic = as.numeric(py_result$waic$pointwise$p_waic_i),
        waic = as.numeric(py_result$waic$pointwise$waic_i)
      )
      
      attr(waic_df, "dims") <- c(n_obs = py_result$waic$n_obs, n_samples = py_result$waic$n_samples)
      class(waic_df) <- c("because_waic", "data.frame")
      
      result$WAIC <- waic_df
    }
    # Compute summary statistics
    sum_stats <- if (!is.null(mcmc_samples)) summary(mcmc_samples) else NULL
    if (!is.null(mcmc_samples) && n.chains > 1) {
      tryCatch({
        n_ch <- length(mcmc_samples)
        first_chain <- as.matrix(mcmc_samples[[1]])
        pnames <- colnames(first_chain)
        n_params <- length(pnames)
        rhat_vals <- numeric(n_params)
        n_iter <- nrow(first_chain)
        for (p in 1:n_params) {
          chain_means <- numeric(n_ch)
          chain_vars <- numeric(n_ch)
          for (c in 1:n_ch) {
            vals <- as.matrix(mcmc_samples[[c]])[, p]
            chain_means[c] <- mean(vals)
            chain_vars[c] <- var(vals)
          }
          grand_mean <- mean(chain_means)
          B <- n_iter * var(chain_means)
          W <- mean(chain_vars)
          if (W > 0) {
            var_plus <- ((n_iter - 1) / n_iter) * W + (1 / n_iter) * B
            rhat_vals[p] <- sqrt(var_plus / W)
          } else {
            rhat_vals[p] <- 1.0
          }
        }
        sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = rhat_vals)
      }, error = function(e) {})
    }
    result$summary <- sum_stats

    result$call <- original_call
    class(result) <- "because"
    return(result)
  }


  # Prune unused tactical variables to avoid JAGS warnings
  for (v in c("N", "zeros", "zero_vec")) {
    if (v %in% names(data) && !grepl(paste0("\\b", v, "\\b"), model_string)) {
      data[[v]] <- NULL
    }
  }

  model_file <- tempfile(fileext = ".jg")
  writeLines(model_string, model_file)

  # If latent variables are present and this is a standard run (dsep=FALSE),
  # print the MAG structure and basis set for user verification, as requested.
  # Display MAG structure for latent variable models (non-dsep runs)
  if (
    !dsep &&
      latent_method == "correlations" &&
      !is.null(latent) &&
      length(latent) > 0 &&
      !quiet
  ) {
    message("--- Latent Variable Structure (MAG) ---")
    if (length(induced_cors) > 0) {
      for (cor_pair in induced_cors) {
        message(sprintf("  Induced Correlation: %s <-> %s", cor_pair[1], cor_pair[2]))
      }
    } else {
      message("  No induced correlations.")
    }
    message("---------------------------------------")
  }

  if (verbose) {
    message("Generated JAGS model:\n", model_string)
  }

  # --- Monitor Parameters ---
  # Handle monitor mode
  monitor_mode <- NULL
  custom_monitors <- character(0)

  if (!is.null(monitor)) {
    if (is.character(monitor)) {
      if ("interpretable" %in% monitor) {
        monitor_mode <- "interpretable"
        custom_monitors <- setdiff(monitor, "interpretable")
      } else if ("all" %in% monitor) {
        monitor_mode <- "all"
        custom_monitors <- setdiff(monitor, "all")
      } else if (identical(monitor, "")) {
        monitor_mode <- "interpretable"
      } else {
        # Entirely custom vector
        custom_monitors <- monitor
      }
    }
  } else {
    monitor_mode <- "interpretable"
  }

  if (
    is.null(monitor) ||
      (!is.null(monitor_mode) && monitor_mode %in% c("interpretable", "all"))
  ) {
    lines <- unlist(strsplit(model_string, "\n"))

    extract_names <- function(pattern) {
      out <- grep(pattern, lines, value = TRUE)
      out <- grep("(<-|~)", out, value = TRUE)
      matches <- regmatches(
        out,
        regexec(
          "(?:logit|log|cloglog|probit)?\\(?\\s*([a-zA-Z0-9_.]+)(?:\\[.*\\])?\\)?\\s*(?:<-|~)",
          out
        )
      )

      res <- sapply(matches, function(m) {
        if (length(m) >= 2) m[2] else NA
      })

      found_params <- as.character(na.omit(res))
      return(found_params)
    }

    # Extract all parameters
    all_params <- unique(c(
      extract_names("^\\s*beta"),
      extract_names("^\\s*alpha"),
      extract_names("^\\s*lambda"),
      extract_names("^\\s*tau"),
      extract_names("^\\s*rho"),
      extract_names("^\\s*sigma"),
      extract_names("^\\s*z_"),
      extract_names("(^|\\W)p_"),
      extract_names("(^|\\W)psi"),
      extract_names("^\\s*r_"),
      extract_names("^\\s*cutpoint")
    ))

    # Ensure mag_exogenous_vars is defined (defaults to empty)
    if (!exists("mag_exogenous_vars")) {
      mag_exogenous_vars <- character(0)
    }

    # Remove tau_obs_* (deterministic constants, not stochastic parameters)
    # Remove tau_obs_* (deterministic constants, not stochastic parameters)
    all_params <- all_params[!grepl("^tau_obs", all_params)]

    if (!is.null(monitor_mode) && monitor_mode == "interpretable") {
      # Filter to interpretable parameters only
      monitor <- all_params[
        (grepl("^alpha", all_params) &
          gsub("^alpha_?", "", all_params) %in% response_vars &
          !gsub("^alpha_?", "", all_params) %in% mag_exogenous_vars) | # Response intercepts, excluding MAG exogenous
          grepl("^beta", all_params) | # All regression coefficients
          grepl("^rho", all_params) | # Induced correlations
          grepl("^sigma", all_params) | # Variance components
          grepl("^psi", all_params) | # Zero-inflation or occupancy probability
          grepl("^p_", all_params) | # Detection probability
          grepl("^z_", all_params) | # Latent state
          grepl("^r_", all_params) | # Negative Binomial size
          grepl("^cutpoint", all_params) | # Ordinal cutpoints
          grepl("^K_", all_params) | # BMA selection indices
          (grepl("^lambda", all_params) &
            gsub("^lambda_?", "", all_params) %in% response_vars &
            !gsub("^lambda_?", "", all_params) %in% mag_exogenous_vars) # Response lambdas, excluding MAG exogenous
      ]
    } else {
      # monitor_mode == "all" or NULL: include everything
      monitor <- all_params

      # Also include response variables (for imputation inspection)
      # Also include response variables (for imputation inspection)
      # We extract them directly from the equations (which include auto-added intercept models)
      response_vars_all <- unique(sapply(equations, function(eq) {
        all.vars(formula(eq)[[2]])
      }))

      # Only add them if they are in the model (obviously)
      if (length(response_vars_all) > 0) {
        # Extension Hook: Map response variables to parameters (e.g. z_Y for occupancy)
        adj_response_vars <- unlist(lapply(response_vars_all, function(v) {
          get_monitor_vars_hook(family_obj, v)
        }))
        monitor <- unique(c(monitor, adj_response_vars))
      }
    }
  }

  # Add custom monitors provided by the user
  if (length(custom_monitors) > 0) {
    if (!is.null(monitor)) {
      monitor <- unique(c(monitor, custom_monitors))
    } else {
      monitor <- custom_monitors
    }
  }

  # Guarantee monitor is not an empty character vector
  if (is.null(monitor) || length(monitor) == 0) {
    if (exists("all_params") && length(all_params) > 0) {
      monitor <- all_params
    }
    if (is.null(monitor) || length(monitor) == 0) {
      if (exists("lines") && length(lines) > 0) {
        assign_lines <- grep("(<-|~)", lines, value = TRUE)
        matches <- regmatches(
          assign_lines,
          regexec("^\\s*([a-zA-Z0-9_.]+)", assign_lines)
        )
        found <- sapply(matches, function(m) if (length(m) >= 2) m[2] else NA)
        monitor <- unique(as.character(na.omit(found)))
      }
    }
  }

  # --- NIMBLE pre-processing ---
  # Ensure all variables used as precision/covariance matrices are numeric matrices
  if (engine == "nimble" && is.list(data)) {
    prec_vars <- names(data)[grepl("^Prec_", names(data))]
    for (pv in prec_vars) {
      if (!is.matrix(data[[pv]])) {
        data[[pv]] <- as.matrix(data[[pv]])
      }
      # [STABILITY] Add small diagonal jitter (nugget) to ensure positive definiteness
      # and prevent numerical singularities during C++ compilation.
      # Ref: Rasmussen & Williams (2006), Gaussian Processes for Machine Learning.
      diag(data[[pv]]) <- diag(data[[pv]]) + 1e-6
    }
  }

  # Add pointwise log-likelihood monitoring if WAIC requested
  # (Future: LOO will also use this)
  if (WAIC) {
    # Extract log_lik parameters from model
    log_lik_params <- unique(c(
      extract_names("^\\s*log_lik")
    ))

    if (length(log_lik_params) > 0) {
      monitor <- unique(c(monitor, log_lik_params))
      if (!quiet) {
        message(
          "Monitoring ",
          length(log_lik_params),
          " pointwise log-likelihood parameter(s) for WAIC"
        )
      }
    }
  }

  # Add response variables
  # Use perl=TRUE for robust regex matching of variable names
  matches <- regmatches(
    model_string,
    gregexpr(
      "\\b([a-zA-Z0-9_.]+)\\s*\\[1:N\\]\\s*~",
      model_string,
      perl = TRUE
    )
  )[[1]]

  response_vars <- unique(gsub("\\s*\\[1:N\\]\\s*~", "", matches))
  for (v in response_vars) {
    # Skip if variable is in variability list (it's latent, not data)
    if (!is.null(variability) && v %in% names(variability_list)) {
      next
    }

    if (!v %in% names(data)) {
      base <- sub("[0-9]+$", "", v)
      if (base %in% names(data)) data[[v]] <- data[[base]]
    }
  }

  # Extension Hook: Custom inits (e.g. occupancy latent states)
  extension_inits <- get_inits_hook(family_obj, data)

  # Clean up data list: Remove variables not present in the model code to avoid warnings
  model_code_str <- model_output$model
  vars_to_remove <- character(0)

  # Always keep these structural/special variables
  keep_vars <- c("zeros")

  for (v in names(data)) {
    if (v %in% keep_vars) {
      next
    }

    # Check if variable appears in the model code as a token
    # Use perl=TRUE for word boundaries
    if (!grepl(paste0("\\b", v, "\\b"), model_code_str, perl = TRUE)) {
      vars_to_remove <- c(vars_to_remove, v)
    }
  }

  if (length(vars_to_remove) > 0) {
    for (v in vars_to_remove) {
      data[[v]] <- NULL
    }
  }

  # Store MCMC compilation and run results
  samples <- NULL
  model <- NULL

  nimble_waic <- NULL
  if (engine == "nimble") {
    # --- NIMBLE EXECUTION PIPELINE ---
    if (!requireNamespace("nimble", quietly = TRUE)) {
      stop(
        "The 'nimble' package is required when engine = 'nimble'.\n",
        "Please install it using: install.packages('nimble')\n",
        "For detailed installation instructions and system requirements (e.g. Rtools/Xcode),\n",
        "see: https://r-nimble.org/download"
      )
    }

    # Attach nimble to the search path to avoid 'getNimbleOption' errors during evaluation
    if (!requireNamespace("nimble", quietly = TRUE)) {
      stop("The 'nimble' package is required for this model but not installed.")
    }

    if (!quiet) {
      message("Compiling model via NIMBLE...")
    }

    # --- NIMBLE Family Optimizations (S3) ---
    # Extensions can implement nimble_family_optimization to provide
    # specialized distributions (e.g. dImperfect) and model transformations.
    nimble_funcs <- list()
    if (!requireNamespace("nimble", quietly = TRUE)) {
      stop("Package 'nimble' is required for engine = 'nimble'.")
    }
    
    if (!("package:nimble" %in% search())) {
      suppressPackageStartupMessages(attachNamespace("nimble"))
    }


    nimble_string <- model_string

    # Generic cleanup for NIMBLE: strip JAGS-specific log-density nodes
    lines <- strsplit(nimble_string, "\n")[[1]]
    lines <- lines[!grepl("logdensity\\.", lines)]
    lines <- lines[!grepl("log_lik_", lines)]
    lines <- lines[!grepl("lik_matrix_", lines)]
    nimble_string <- paste(lines, collapse = "\n")

    for (v in names(family)) {
      fam_obj <- structure(
        list(name = family[[v]]),
        class = c(paste0("because_family_", family[[v]]), "because_family")
      )
      opt_res <- nimble_family_optimization(
        fam_obj,
        nimble_string,
        variable = v
      )
      nimble_string <- opt_res$model_string
      if (length(opt_res$nimble_functions) > 0) {
        nimble_funcs <- c(nimble_funcs, opt_res$nimble_functions)
      }

      # If discretized latent states were marginalized, remove them from monitors
      if (
        nimble_string != model_string && any(grepl(paste0("z_", v), monitor))
      ) {
        monitor <- setdiff(monitor, paste0("z_", v))
      }
    }

    # Register nimble functions to local environment for compiler
    if (length(nimble_funcs) > 0) {
      unique_names <- unique(names(nimble_funcs))
      for (fn_name in unique_names) {
        # Assign to local environment so nimbleModel can find them
        assign(fn_name, nimble_funcs[[fn_name]], envir = environment())

        # Explicitly register with NIMBLE if it looks like a distribution
        if (startsWith(fn_name, "d")) {
          try(
            nimble::registerDistributions(nimble_funcs[fn_name]),
            silent = TRUE
          )
        }
      }
    }

    # Ensure all monitored parameters have initial values for NIMBLE stability
    # JAGS auto-initializes many nodes, but NIMBLE is more rigorous.
    # Latent random effects and categorical intercepts must be initialized.
    nimble_inits <- extension_inits
    if (is.null(nimble_inits)) nimble_inits <- list()
    for (p in monitor) {
      if (!p %in% names(nimble_inits)) {
        if (grepl("^(tau_|sigmay_|sigmap_|sigmar_|sigma_)", p)) {
          nimble_inits[[p]] <- 1.0 # Standard unit variance start
        } else if (grepl("^alpha_", p)) {
          # [STABILITY] Initialize intercepts closer to data mean if available
          resp_name <- sub("^alpha_", "", p)
          if (resp_name %in% names(data)) {
              m_val <- mean(as.numeric(data[[resp_name]]), na.rm = TRUE)
              # If it looks like a count or binary, use link Scale
              if (all(as.numeric(data[[resp_name]]) >= 0, na.rm = TRUE)) {
                  nimble_inits[[p]] <- log(max(0.1, m_val))
              } else {
                  nimble_inits[[p]] <- m_val
              }
          } else {
            nimble_inits[[p]] <- 0.0
          }
        } else if (grepl("^beta_", p)) {
          # [STABILITY] Start slopes at 0 but ensure they will be jittered
          nimble_inits[[p]] <- 0.0
        } else if (grepl("^psi_", p)) {
          nimble_inits[[p]] <- 0.5
        } else if (grepl("^r_", p)) {
          nimble_inits[[p]] <- 1.0
        } else if (grepl("^sigma_total_", p)) {
          # [PARTITIONING] sigma_total drives both tau_u and tau_res via lambda.
          # Starting at 0 collapses the dmnorm prior to a point mass (tau -> Inf).
          nimble_inits[[p]] <- 1.0
        } else if (grepl("^lambda_", p)) {
          nimble_inits[[p]] <- 0.5
        } else if (grepl("^cutpoint", p)) {
          # Ordinal cutpoints need to be ordered; leaving as 0 can crash
          # better to use extension_inits or stay conservative
          nimble_inits[[p]] <- 0.0
        }
      }
    }

    # Convert the JAGS model string directly into a NIMBLE model
    nimble_model <- tryCatch(
      {
        nimble_string <- sub("^\\s*model\\s*\\{", "{", nimble_string)
        nimble_code <- parse(text = nimble_string)[[1]]

        nimble_constants <- data
        nimble_data <- list()
        if (!is.null(data[["L_multiPhylo"]])) {
            nimble_data[["L_multiPhylo"]] <- data[["L_multiPhylo"]]
            nimble_constants[["L_multiPhylo"]] <- NULL
        }
        if (!is.null(data[["Prec_multiPhylo"]])) {
            nimble_data[["Prec_multiPhylo"]] <- data[["Prec_multiPhylo"]]
            nimble_constants[["Prec_multiPhylo"]] <- NULL
        }
        if (!is.null(data[["L_phylo"]])) {
            # Single-tree non-centered Cholesky factor: treat as data (observed matrix)
            nimble_data[["L_phylo"]] <- data[["L_phylo"]]
            nimble_constants[["L_phylo"]] <- NULL
        }

        m_obj <- suppressMessages(suppressWarnings(nimble::nimbleModel(
          code = nimble_code,
          constants = nimble_constants,
          data = nimble_data,
          inits = nimble_inits,
          buildDerivs = (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC"),
          calculate = FALSE
        )))
        
        # [NEW] Ensure all parameters are initialized properly for NIMBLE
        # This handles vector/matrix nodes like alpha_y and err_y
        model_nodes <- m_obj$getNodeNames(stochOnly = TRUE, includeData = FALSE)
        for (node in model_nodes) {
          # Only initialize if NOT already set in nimble_inits
          # (getNodeNames returns specific indices, so we strip them)
          base_node <- sub("\\[.*\\]", "", node)
          if (!any(grepl(paste0("^", base_node, "$"), names(nimble_inits)))) {
            # Check if this is a variance/precision node
            is_prec <- grepl("^(tau_|sigmay_|sigmap_|sigmar_)", node)
            val <- if (is_prec) 1 else 0
            if (grepl("^psi_", node)) val <- 0.5
            if (grepl("^r_", node)) val <- 1
            # [PARTITIONING] sigma_total_ and lambda_ MUST NOT start at 0.
            # sigma_total=0 -> tau=Inf -> dmnorm at point mass -> NA cascade.
            if (grepl("^sigma_total_", node)) val <- 1.0
            if (grepl("^lambda_", node)) val <- 0.5
            
            try({
              curr_val <- m_obj[[node]]
              if (any(is.na(curr_val)) || any(is.nan(curr_val))) {
                m_obj[[node]] <- val
              }
            }, silent = TRUE)
          }
        }
        
        # [CRITICAL FIX] Export the auto-initialized values back into nimble_inits
        # so they can be sent to parallel workers.
        unique_base_nodes <- unique(sub("\\[.*\\]", "", model_nodes))
        for (v in unique_base_nodes) {
            if (!v %in% names(nimble_inits)) {
                try({
                    nimble_inits[[v]] <- m_obj[[v]]
                }, silent = TRUE)
            }
        }
        
        m_obj
      },
      error = function(e) {
        if (!quiet) {
          message("\nCRITICAL NIMBLE ERROR during model initialization:")
          message(e$message)
        }
        stop(paste(e, "\n\n", model_string))
      }
    )

    # Configure MCMC
    mcmc_conf <- nimble::configureMCMC(
      nimble_model,
      monitors = monitor,
      enableWAIC = WAIC
    )

    # Check if HMC was specifically requested
    if (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC") {
      if (!requireNamespace("nimbleHMC", quietly = TRUE)) {
        stop("The 'nimbleHMC' package is required to use HMC samplers in NIMBLE. Install it with install.packages('nimbleHMC')")
      }
      if (!("package:nimbleHMC" %in% search())) {
        suppressPackageStartupMessages(attachNamespace("nimbleHMC"))
      }
      if (!quiet) {
        message("Building and compiling NIMBLE MCMC using nimbleHMC::buildHMC (this may take a moment)...")
      }
      nimble_mcmc <- nimbleHMC::buildHMC(nimble_model)
    } else {
      # Harden NIMBLE sampler assignments.
      because::nimble_harden_samplers(
        mcmc_conf,
        family          = family,
        nimble_samplers = nimble_samplers,
        quiet           = quiet
      )

      if (!quiet) {
        message("Building and compiling NIMBLE MCMC (this may take a moment)...")
      }
      nimble_mcmc <- nimble::buildMCMC(mcmc_conf)
    }
    
    # [PERFORMANCE] Skip main-thread compilation if running in parallel
    # to avoid redundant triple/quadruple compilation overhead.
    if (!parallel || n.cores == 1 || n.chains == 1) {
      compiled_model <- nimble::compileNimble(nimble_model)
      compiled_mcmc <- nimble::compileNimble(nimble_mcmc, project = nimble_model)
    } else {
      compiled_mcmc <- NULL
    }

    if (parallel && n.cores > 1 && n.chains > 1) {
      # Parallel execution for NIMBLE
      if (!quiet) {
        message(sprintf(
          "Running %d NIMBLE chains in %s on %d cores...",
          n.chains,
          "parallel",
          n.cores
        ))
      }

      if (is.null(cl)) {
        cl <- parallel::makeCluster(n.cores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      }

      # Helper for parallel NIMBLE chain
      run_nimble_chain <- function(
        chain_id,
        model_string,
        data,
        family,
        nimble_inits, # Corrected: Pass populated inits instead of extension_inits
        monitor,
        n.iter,
        n.burnin,
        n.thin,
        WAIC,
        nimble_samplers,
        quiet
      ) {
        if (!requireNamespace("nimble", quietly = TRUE)) {
          return(NULL)
        }
        if (!requireNamespace("nimble", quietly = TRUE)) {
          stop(
            "The 'nimble' package is required for this model but not installed."
          )
        }

        # --- NIMBLE Family Optimizations via S3 (Worker) ---
        nimble_funcs <- list()
        nimble_string <- model_string

        # Cleanup JAGS nodes
        lines <- strsplit(nimble_string, "\n")[[1]]
        lines <- lines[!grepl("logdensity\\.", lines)]
        lines <- lines[!grepl("log_lik_", lines)]
        lines <- lines[!grepl("lik_matrix_", lines)]
        nimble_string <- paste(lines, collapse = "\n")

        for (v in names(family)) {
          fam_obj <- structure(
            list(name = family[[v]]),
            class = c(paste0("because_family_", family[[v]]), "because_family")
          )
          opt_res <- nimble_family_optimization(
            fam_obj,
            nimble_string,
            variable = v
          )
          nimble_string <- opt_res$model_string
          if (length(opt_res$nimble_functions) > 0) {
            nimble_funcs <- c(nimble_funcs, opt_res$nimble_functions)
          }
          if (
            nimble_string != model_string &&
              any(grepl(paste0("z_", v), monitor))
          ) {
            monitor <- setdiff(monitor, paste0("z_", v))
          }
        }
        # Assign functions to local environment so nimbleModel can find them
        if (length(nimble_funcs) > 0) {
          unique_names <- unique(names(nimble_funcs))
          for (fn_name in unique_names) {
            assign(fn_name, nimble_funcs[[fn_name]], envir = environment())
            
            # Explicitly register with NIMBLE if it looks like a distribution
            if (startsWith(fn_name, "d")) {
              try(nimble::registerDistributions(nimble_funcs[fn_name]), silent = TRUE)
            }
          }
        }

        # Strip model { ... } wrapping
        nimble_string <- sub("^\\s*model\\s*\\{", "{", nimble_string)
        nimble_code <- parse(text = nimble_string)[[1]]

        # [STABILITY] Jitter inits for this specific chain
        curr_inits <- nimble_inits # Corrected: Use populated inits
        if (is.null(curr_inits)) curr_inits <- list()
        
        # Add a stochastic jitter to all continuous parameters
        # This is critical for NIMBLE to escape locally flat regions
        # We increase the range to 0.1 for more robust exploration
        set.seed(12345 + chain_id)
        for (p_name in names(curr_inits)) {
            val <- curr_inits[[p_name]]
            if (is.numeric(val) && length(val) == 1) {
                if (grepl("^(beta_|alpha_)", p_name)) {
                    curr_inits[[p_name]] <- val + rnorm(1, 0, 0.1)
                } else if (grepl("^(tau_|sigma_)", p_name)) {
                    curr_inits[[p_name]] <- max(0.1, val * exp(rnorm(1, 0, 0.1)))
                }
            }
        }

        # Ensure NIMBLE namespace is loaded on the worker
        requireNamespace("nimble", quietly = TRUE)
        if (!("package:nimble" %in% search())) {
          suppressPackageStartupMessages(attachNamespace("nimble"))
        }
        
        # If HMC is requested, ensure nimbleHMC is attached so samplers are registered
        if (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC") {
          requireNamespace("nimbleHMC", quietly = TRUE)
          if (!("package:nimbleHMC" %in% search())) {
            suppressPackageStartupMessages(attachNamespace("nimbleHMC"))
          }
        }

        nimble_constants <- data
        nimble_data <- list()
        if (!is.null(data[["L_multiPhylo"]])) {
            nimble_data[["L_multiPhylo"]] <- data[["L_multiPhylo"]]
            nimble_constants[["L_multiPhylo"]] <- NULL
        }
        if (!is.null(data[["Prec_multiPhylo"]])) {
            nimble_data[["Prec_multiPhylo"]] <- data[["Prec_multiPhylo"]]
            nimble_constants[["Prec_multiPhylo"]] <- NULL
        }
        if (!is.null(data[["L_phylo"]])) {
            # Single-tree non-centered Cholesky factor: treat as data (observed matrix)
            nimble_data[["L_phylo"]] <- data[["L_phylo"]]
            nimble_constants[["L_phylo"]] <- NULL
        }

        worker_model <- suppressMessages(suppressWarnings(nimble::nimbleModel(
          code = nimble_code,
          constants = nimble_constants,
          data = nimble_data,
          inits = curr_inits,
          buildDerivs = (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC"),
          calculate = FALSE
        )))

        worker_conf <- nimble::configureMCMC(
          worker_model,
          monitors = monitor,
          enableWAIC = WAIC
        )
        
        # Harden worker MCMC samplers via the installed package function.
        if (is.character(nimble_samplers) && length(nimble_samplers) == 1 && nimble_samplers == "HMC") {
          if (!requireNamespace("nimbleHMC", quietly = TRUE)) {
            stop("The 'nimbleHMC' package is required to use HMC samplers in NIMBLE.")
          }
          worker_mcmc <- nimbleHMC::buildHMC(worker_model)
        } else {
          because::nimble_harden_samplers(
            worker_conf,
            family          = family,
            nimble_samplers = nimble_samplers,
            quiet           = TRUE
          )
          worker_mcmc <- nimble::buildMCMC(worker_conf)
        }
        worker_c_model <- nimble::compileNimble(worker_model)
        worker_c_mcmc <- nimble::compileNimble(
          worker_mcmc,
          project = worker_model
        )

        # --- Execute MCMC ---
        # Wrap everything in try to avoid hanging the cluster on error
        res <- try({
          samples <- nimble::runMCMC(
            worker_c_mcmc,
            niter = n.iter,
            nburnin = n.burnin,
            nchains = 1,
            thin = n.thin,
            samplesAsCodaMCMC = TRUE,
            WAIC = WAIC
          )
          samples
        }, silent = TRUE)

        if (inherits(res, "try-error")) {
          return(paste("NIMBLE WORKER ERROR:", as.character(res)))
        }
        return(res)
      }

      # Export all required variables to worker nodes
      parallel::clusterExport(
        cl,
        c(
          "model_string", "data", "family", "nimble_inits",
          "monitor", "n.iter", "n.burnin", "n.thin",
          "WAIC", "quiet", "run_nimble_chain", "nimble_samplers"
        ),
        envir = environment()
      )

      # Copy S3 methods to workers if package isn't installed
      # (Necessary for devtools::load_all sessions)
      if ("package:because" %in% search()) {
        parallel::clusterEvalQ(cl, library(because))
      }
      if ("package:because.phybase" %in% search()) {
        parallel::clusterEvalQ(cl, library(because.phybase))
      }

      chain_results <- parallel::parLapply(cl, seq_len(n.chains), function(i) {
        run_nimble_chain(
          chain_id = i,
          model_string = model_string,
          data = data,
          family = if (!is.null(family)) as.list(family) else NULL,
          nimble_inits = nimble_inits,
          monitor = monitor,
          n.iter = n.iter,
          n.burnin = n.burnin,
          n.thin = n.thin,
          WAIC = WAIC,
          nimble_samplers = nimble_samplers,
          quiet = quiet
        )
      })

      # Check for worker errors
      for (i in seq_along(chain_results)) {
        if (is.character(chain_results[[i]])) {
          stop(sprintf("Chain %d failed: %s", i, chain_results[[i]]))
        }
      }

      # Extract NIMBLE WAIC if available
      if (WAIC && is.list(chain_results[[1]]) && !is.null(chain_results[[1]]$WAIC)) {
          nimble_waic <- chain_results[[1]]$WAIC
      }

      # Format into mcmc.list
      samples <- coda::mcmc.list(lapply(chain_results, function(x) {
          if (is.list(x) && !is.null(x$samples)) x$samples else x
      }))
    } else {
      # Sequential execution (Existing logic)
      if (!quiet) {
        message(sprintf(
          "Sampling %d chains sequentially via NIMBLE...",
          n.chains
        ))
      }

      # --- Sequential NIMBLE execution via run_nimble_model() ---
      nimble_run_res <- run_nimble_model(
        nimble_code  = nimble_code,
        constants    = data,
        data         = list(),
        inits        = nimble_inits,
        n.chains     = n.chains,
        n.iter       = n.iter,
        n.burnin     = n.burnin,
        n.thin       = n.thin,
        monitor      = monitor,
        quiet        = quiet,
        family       = family,
        nimble_samplers = nimble_samplers,
        nimble_waic  = WAIC
      )
      samples     <- nimble_run_res$samples
      nimble_waic <- nimble_run_res$WAIC
    }
    model <- nimble_model # Store the R (uncompiled) model object
    # Store compiled objects for because_continue() --- NULL when parallel=TRUE
    saved_nimble_compiled <- if (exists("compiled_mcmc")) compiled_mcmc else NULL
    saved_nimble_cmodel   <- if (exists("compiled_model")) compiled_model else NULL
    saved_nimble_samplers <- nimble_samplers
  } else {
    # --- JAGS EXECUTION PIPELINE (Default) ---
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
        # Load rjags in each worker
        if (!requireNamespace("rjags", quietly = TRUE)) {
          stop("Package 'rjags' is required for parallel execution.")
        }
        # Explicitly load rjags to ensure modules are available
        loadNamespace("rjags")

        # Compile model for this chain
        # Explicitly set RNG seed to ensure chains are different
        inits_list <- c(
          extension_inits,
          list(
            .RNG.name = "base::Wichmann-Hill",
            .RNG.seed = 12345 + chain_id
          )
        )

        par_inits <- list(inits_list)
        model <- run_jags_model(model_file, data, par_inits, 1L, n.adapt, quiet, model_string)

        # Burn-in
        if (n.burnin > 0) {
          update(model, n.iter = n.burnin)
        }

        samples <- sample_jags_model(model, monitor, n.iter, n.burnin, n.thin)

        return(list(samples = samples, model = model))
      }

      # Export necessary objects to cluster, including backend helpers
      parallel::clusterExport(cl, c("run_single_chain", "run_jags_model", "sample_jags_model"), envir = environment())
      # Ensure because package (and its helpers) are available on workers
      parallel::clusterEvalQ(cl, {
        if (requireNamespace("because", quietly = TRUE)) library(because)
      })

      # Run chains in parallel
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

      # Combine samples from all chains
      if (!is.null(chain_results[[1]]$samples)) {
        samples <- coda::mcmc.list(lapply(chain_results, function(x) {
          x$samples[[1]]
        }))
      } else {
        samples <- NULL
      }

      # Use the first chain's model for DIC/WAIC (they all have the same structure)
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

      # Combine with occupancy inits
      inits_list <- lapply(1:n.chains, function(i) {
        c(
          extension_inits,
          list(
            .RNG.name = "base::Wichmann-Hill",
            .RNG.seed = 12345 + i
          )
        )
      })

      # Model file is managed by tryCatch below (diagostic dump on failure)

      model <- run_jags_model(model_file, data, inits_list, n.chains, n.adapt, quiet, model_string)
      if (n.burnin > 0) {
        update(model, n.iter = n.burnin)
      }

      # Disable DIC/WAIC if only 1 chain (rjags requirement)
      if (n.chains < 2 && (DIC || WAIC)) {
        warning(
          "DIC and WAIC require at least 2 chains. Disabling calculation."
        )
        DIC <- FALSE
        WAIC <- FALSE
      }

      samples <- sample_jags_model(model, monitor, n.iter, n.burnin, n.thin)
    }
  }

  # Summarize posterior
  if (!is.null(samples)) {
    sum_stats <- summary(samples)
  } else {
    sum_stats <- NULL
  }

  # Explicitly calculate R-hat if multiple chains
  if (!is.null(samples) && n.chains > 1) {
    tryCatch(
      {
        # Manual R-hat calculation to avoid coda::gelman.diag issues
        # with parallel chains and R scoping problems
        n_chains <- length(samples)

        # Use base R colnames to get parameter names
        first_chain <- as.matrix(samples[[1]])
        pnames <- colnames(first_chain)
        n_params <- length(pnames)

        psrf <- matrix(NA, nrow = n_params, ncol = 2)
        rownames(psrf) <- pnames
        colnames(psrf) <- c("Point est.", "Upper C.I.")

        # Convert chains to matrices ONCE outside the loop
        chain_matrices <- lapply(samples, as.matrix)

        for (idx in seq_len(n_params)) {
          # Extract column idx from each chain matrix
          vals <- do.call(cbind, lapply(chain_matrices, function(m) m[, idx]))

          # Check for constant chains (variance 0)
          # Use explicit variance calculation to avoid R scoping issues
          chain_vars <- numeric(ncol(vals))
          for (col_idx in seq_len(ncol(vals))) {
            chain_vars[col_idx] <- stats::var(vals[, col_idx])
          }

          if (any(chain_vars < 1e-10)) {
            psrf[idx, 1] <- 1.0 # If constant, Rhat is 1
            next
          }

          # Calculate B/W (Gelman-Rubin statistic)
          n_samples <- nrow(vals)
          chain_means <- colMeans(vals)

          # Between-chain variance
          B <- n_samples * stats::var(chain_means)

          # Within-chain variance
          W <- mean(chain_vars)

          # Estimated variance
          var_plus <- (n_samples - 1) / n_samples * W + B / n_samples

          # R-hat
          rhat <- sqrt(var_plus / W)
          psrf[idx, 1] <- rhat
        }
        # Add R-hat to summary statistics
        # summary(samples) returns a list with 'statistics' and 'quantiles'
        # We want to add R-hat to the statistics matrix

        # Match parameter names
        common_params <- intersect(
          rownames(sum_stats$statistics),
          rownames(psrf)
        )

        if (length(common_params) > 0) {
          sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = NA)
          rhat_col_idx <- which(colnames(sum_stats$statistics) == "Rhat")

          # Use numeric indexing to avoid any strange symbol evaluation
          for (j in seq_along(common_params)) {
            p <- common_params[j]
            row_idx <- which(rownames(sum_stats$statistics) == p)
            psrf_row_idx <- which(rownames(psrf) == p)
            if (length(row_idx) == 1 && length(psrf_row_idx) == 1) {
              sum_stats$statistics[row_idx, rhat_col_idx] <- psrf[
                psrf_row_idx,
                1
              ]
            }
          }
        }
      },
      error = function(e) {
        warning("Could not calculate R-hat: ", e$message)
      }
    )
  }

  # Filter internal parameters (log_lik) from summary parameters
  # We keep them in samples for WAIC calculation but hide them from the summary output
  if (!is.null(sum_stats)) {
    if (is.matrix(sum_stats$statistics)) {
      rows_to_keep <- !grepl("^log_lik", rownames(sum_stats$statistics))
      sum_stats$statistics <- sum_stats$statistics[
        rows_to_keep,
        ,
        drop = FALSE
      ]
      sum_stats$quantiles <- sum_stats$quantiles[rows_to_keep, , drop = FALSE]
    } else {
      # Single parameter case (statistics is a vector)
      # Check if the single parameter is log_lik
      param_name <- colnames(samples[[1]])
      if (length(param_name) == 1 && grepl("^log_lik", param_name)) {
        # If the only parameter is log_lik, return empty stats
        # Or handle appropriately. For now, empty seems safest or just nullify.
        sum_stats <- NULL
      }
    }
  }

  # Initialize result object
  result <- list(
    model = model,
    model_code = model_output$model,
    data = data, # Store data for recompilation if needed
    input = list(
      equations = equations,
      random = random,
      structure = structure,
      data = original_data, # Store original data too for safety
      latent = latent,
      distribution = distribution,
      family = if (!is.null(family)) as.list(family) else NULL,
      variability = variability,
      poly_terms = all_poly_terms # Needed by plot_dag to reconstruct diamond nodes
    ),
    samples = samples,
    summary = sum_stats,
    monitor = monitor,
    modfile = model_file,
    dsep = dsep,
    dsep_tests = dsep_tests,
    dsep_results = dsep_results,
    parameter_map = parameter_map,
    induced_correlations = induced_cors,
    scale_info = scale_info,
    stacked_data = if (exists("stack_res") && stack_res$is_stacked) {
      stack_res$data
    } else {
      NULL
    }
  )

  # Combine with basic model info
  # NOTE: result$model already holds the live rjags object (from result list above).
  # result$model_code holds the model string for reference / recompilation.
  # Do NOT overwrite result$model here --- keeping the live object enables because_continue().
  result$engine         <- engine
  if (engine == "nimble") {
    result$nimble_compiled <- if (exists("saved_nimble_compiled")) saved_nimble_compiled else NULL
    result$nimble_cmodel   <- if (exists("saved_nimble_cmodel")) saved_nimble_cmodel else NULL
    result$nimble_samplers <- if (exists("saved_nimble_samplers")) saved_nimble_samplers else NULL
  }
  result$samples        <- samples
  result$parameter_map  <- parameter_map
  result$data           <- data
  result$original_data  <- original_data
  result$family         <- family
  result$categorical_vars <- attr(data, "categorical_vars")
  result$poly_terms     <- all_poly_terms
  result$equations      <- equations
  result$parallel       <- parallel
  result$n.cores        <- n.cores

  # --- Result Enrichment ---
  # If we have a structure, try to extract labels for ordering
  if (!is.null(structure)) {
    result$species_order <- get_order_labels_hook(structure)
  } else if (!is.null(id_col) && is.data.frame(original_data)) {
    # If no structure but ID col provided
    result$species_order <- as.character(original_data[[id_col]])
  }

  # Assign class immediately (needed for print/summary/waic methods)
  result$call <- original_call
  class(result) <- "because"

  # Add DIC and WAIC
  # For parallel runs, recompile the model if ic_recompile=TRUE
  if (
    (DIC || WAIC) && parallel && n.cores > 1 && n.chains > 1 && ic_recompile
  ) {
    message("Recompiling model for DIC/WAIC calculation...")

    # Recompile model with 2 chains for IC calculation (DIC requires >=2)
    ic_inits <- lapply(1:2, function(i) {
      c(
        extension_inits,
        list(
          .RNG.name = "base::Wichmann-Hill",
          .RNG.seed = 12345 + i
        )
      )
    })

    ic_model <- rjags::jags.model(
      model_file,
      data = data,
      inits = ic_inits,
      n.chains = 2,
      n.adapt = n.adapt,
      quiet = quiet
    )

    # Short burn-in (use a fraction of original)
    if (n.burnin > 0) {
      update(ic_model, n.iter = min(n.burnin, 500))
    }

    # Compute DIC
    if (DIC) {
      if (n.iter > n.burnin) {
        result$DIC <- rjags::dic.samples(
          ic_model,
          n.iter = min(n.iter - n.burnin, 1000)
        )
      } else {
        result$DIC <- NULL
      }
    }

    # Compute WAIC using pointwise log-likelihoods
    # Note: WAIC calculation generally uses the posterior samples already collected.
    # We defer WAIC calculation to the common block at the end to ensure consistency.
    # if (WAIC) {
    #   result$WAIC <- because_waic(result)
    # }
  } else if ((DIC || WAIC) && parallel && n.cores > 1 && n.chains > 1) {
    # Parallel without recompilation - warn user
    if (DIC) {
      warning(
        "DIC calculation disabled for parallel chains. Set ic_recompile=TRUE to compute DIC."
      )
      result$DIC <- NULL
    }
    if (WAIC) {
      # WAIC can be computed from pointwise log-likelihoods even with parallel chains
      # Defer to common block
      # result$WAIC <- because_waic(result)
    }
  } else {
    # Sequential execution - use standard approach
  if (DIC) {
    if (engine == "jags") {
        if (n.iter > n.burnin) {
          result$DIC <- rjags::dic.samples(model, n.iter = n.iter - n.burnin)
        } else {
          result$DIC <- NULL
        }
      } else {
        result$DIC <- NULL # NIMBLE does not use rjags::dic.samples
      }
    }
    # WAIC will be computed after class assignment
  }

  # Assign class before WAIC computation (because_waic needs this)
  # Already assigned earlier
  # class(result) <- "because"

  # Compute WAIC if requested (must be after class assignment)
  if (WAIC) {
    if (engine == "nimble" && !is.null(nimble_waic)) {
       # NIMBLE built-in WAIC: lppd = log pointwise predictive density (no penalty),
       # pWAIC = effective parameters, WAIC = -2*(lppd - pWAIC).
       # elpd_waic (as in because) = lppd - pWAIC  (NOT just lppd)
       nimble_elpd     <- nimble_waic$lppd - nimble_waic$pWAIC
       nimble_pwaic    <- nimble_waic$pWAIC
       nimble_waic_val <- nimble_waic$WAIC   # = -2 * nimble_elpd
       # n_obs: N observations x number of modelled response variables
       n_obs_nimble <- tryCatch({
           n_resp <- length(unique(result$parameter_map$response))
           as.integer(data$N * n_resp)
       }, error = function(e) NA_integer_)
       n_samples_nimble <- as.integer((n.iter - n.burnin) / n.thin) * n.chains
       waic_df <- data.frame(
           Estimate = c(nimble_elpd, nimble_pwaic, nimble_waic_val),
           SE = c(NA_real_, NA_real_, NA_real_),  # no pointwise SE from NIMBLE built-in WAIC
           row.names = c("elpd_waic", "p_waic", "waic")
       )
       attr(waic_df, "dims") <- c(n_obs = n_obs_nimble, n_samples = n_samples_nimble)
       class(waic_df) <- c("because_waic", "data.frame")
       result$WAIC <- waic_df
    } else {
       result$WAIC <- because_waic(result)
    }
  }

  # Preserve hierarchical metadata for diagnostics even if data was flat
  result$hierarchical_info <- hierarchical_info
  
  return(result)
}

