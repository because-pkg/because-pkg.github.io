

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
    because_py <- setup_numpyro_environment(parallel = parallel, n.cores = n.cores, n.chains = n.chains)
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

  # --- Selective Variable Collection & Categorical Preprocessing ---
  filtered_res <- collect_and_filter_model_data(
    data           = data,
    equations      = equations,
    random         = random,
    id_col         = id_col,
    link_vars      = link_vars,
    family         = family,
    variability    = variability,
    expand_ordered = expand_ordered,
    quiet          = quiet
  )
  data                                <- filtered_res$data
  original_raw_data_before_preprocess <- filtered_res$original_raw_data_before_preprocess
  fixed_eqs_temp                      <- filtered_res$fixed_eqs_temp
  parsed_random_temp                  <- filtered_res$parsed_random_temp
  global_random_vars                  <- filtered_res$global_random_vars

  # --- Hierarchical Data Detection & Validation ---
  hier_res <- detect_and_validate_hierarchical(
    data               = data,
    equations          = equations,
    fixed_eqs_temp     = fixed_eqs_temp,
    parsed_random_temp = parsed_random_temp,
    global_random_vars = global_random_vars,
    levels             = levels,
    multiscale         = multiscale,
    hierarchy          = hierarchy,
    link_vars          = link_vars,
    latent             = latent,
    random             = random,
    structure          = structure,
    structure_multi    = structure_multi,
    structure_levels   = structure_levels,
    id_col             = id_col,
    quiet              = quiet
  )
  is_hierarchical   <- hier_res$is_hierarchical
  hierarchical_info <- hier_res$hierarchical_info
  levels            <- hier_res$levels
  multiscale        <- hier_res$multiscale
  link_vars         <- hier_res$link_vars
  row_ids           <- hier_res$row_ids

  # --- Random Effects Parsing ---
  rand_res     <- parse_model_random_effects(equations, random)
  equations    <- rand_res$equations
  random_terms <- rand_res$random_terms

  # --- Polynomial Term Extraction ---
  poly_res <- process_polynomial_terms(
    equations         = equations,
    is_hierarchical   = is_hierarchical,
    levels            = levels,
    hierarchical_info = hierarchical_info,
    quiet             = quiet
  )
  equations         <- poly_res$equations
  all_poly_terms    <- poly_res$all_poly_terms
  levels            <- poly_res$levels
  hierarchical_info <- poly_res$hierarchical_info

  # Initialize result variables
  dsep_tests        <- NULL
  dsep_results      <- NULL
  dsep_correlations <- NULL
  induced_cors      <- NULL

  # Handle global variability setting (e.g. variability = "reps")
  variability <- normalize_global_variability(variability, equations, random, id_col)

  # --- Data Frame Preprocessing ---
  original_data <- data

  # --- Hierarchical Data Assembly ---
  if (is_hierarchical) {
    assem_res <- assemble_hierarchical_dataset(
      data                                = data,
      hierarchical_info                   = hierarchical_info,
      equations                           = equations,
      random_terms                        = random_terms,
      latent                              = latent,
      all_poly_terms                      = all_poly_terms,
      original_raw_data_before_preprocess = original_raw_data_before_preprocess,
      engine                              = engine,
      quiet                               = quiet
    )
    data              <- assem_res$data
    hierarchical_info <- assem_res$hierarchical_info
    original_data     <- assem_res$original_data
  }

  # --- Data Validation ---
  validate_binomial_responses(data, equations, family)

  # --- Random Effects Data Prep ---
  re_res <- prepare_random_effects_data(
    data = data, random_terms = random_terms, equations = equations,
    hierarchical_info = hierarchical_info, is_hierarchical = is_hierarchical,
    levels = levels, family = family, quiet = quiet,
    variability = variability, id_col = id_col,
    all_poly_terms = all_poly_terms, latent = latent, structure = structure,
    original_data = original_data
  )
  data              <- re_res$data
  random_structures <- re_res$random_structures
  if (!is.null(re_res$hierarchical_info)) hierarchical_info <- re_res$hierarchical_info
  if (!is.null(re_res$original_data)) original_data <- re_res$original_data

  # --- Structure Processing ---
  struct_res <- process_because_structures(
    structure = structure, structures = structures, structure_obj = structure_obj,
    data = data, hierarchical_info = hierarchical_info,
    is_hierarchical = is_hierarchical, equations = equations,
    family = family, family_obj = family_obj, levels = levels,
    multiscale = multiscale, latent = latent, quiet = quiet,
    engine = engine,
    original_raw_data_before_preprocess = original_raw_data_before_preprocess,
    row_ids = row_ids
  )
  data              <- struct_res$data
  structures        <- struct_res$structures
  is_multiple       <- struct_res$is_multiple
  if (!is.null(struct_res$hierarchical_info)) hierarchical_info <- struct_res$hierarchical_info
  N                 <- struct_res$N

  # --- Missing Data, Variability, and Latent Processing ---
  mv_res <- process_missing_and_variability(
    data                  = data,
    hierarchical_info     = hierarchical_info,
    N                     = N,
    structures            = structures,
    family                = family,
    family_obj            = family_obj,
    equations             = equations,
    fix_residual_variance = fix_residual_variance,
    variability           = variability,
    latent                = latent,
    all_poly_terms        = all_poly_terms,
    quiet                 = quiet,
    dsep                  = dsep,
    latent_method         = latent_method
  )
  data                  <- mv_res$data
  family                <- mv_res$family
  fix_residual_variance <- mv_res$fix_residual_variance
  equations             <- mv_res$equations
  variability           <- mv_res$variability
  variability_list      <- mv_res$variability_list
  latent                <- mv_res$latent
  response_vars_with_na <- mv_res$response_vars_with_na


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
    prior_scale_fixed = prior_scale_fixed, verbose = verbose,
    family_obj = family_obj, hierarchy = hierarchy,
    original_call = original_call, original_data = original_data,
    random = random, response_vars_with_na = response_vars_with_na,
    reuse_models = reuse_models, structure = structure,
    because_py = if (exists("because_py")) because_py else NULL,
    WAIC = WAIC, DIC = DIC
  )
  if (inherits(dsep_res, "because")) {
    return(dsep_res)
  }
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
    return(run_numpyro_pipeline(
      equations = equations,
      data = data,
      family = family,
      structures = structures,
      priors = priors,
      random_terms = random_terms,
      hierarchical_info = hierarchical_info,
      is_hierarchical = is_hierarchical,
      latent = latent,
      n.chains = n.chains,
      n.iter = n.iter,
      n.burnin = n.burnin,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      prior_scale_fixed = prior_scale_fixed,
      quiet = quiet,
      WAIC = WAIC,
      model_string = model_string,
      parameter_map = parameter_map,
      original_call = original_call,
      engine = engine
    ))
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

  # --- Monitor and Inits Preparation ---
  response_vars <- unique(vapply(
    equations,
    function(eq) as.character(all.vars(eq[[2]])[1]),
    character(1)
  ))
  mon_inits <- prepare_monitors_and_inits(
    model_string        = model_string,
    monitor             = monitor,
    response_vars       = response_vars,
    mag_exogenous_vars  = if (exists("mag_exogenous_vars")) mag_exogenous_vars else character(0),
    family_obj          = family_obj,
    equations           = equations,
    engine              = engine,
    data                = data,
    variability         = variability,
    variability_list    = variability_list,
    WAIC                = WAIC,
    quiet               = quiet
  )
  monitor         <- mon_inits$monitor
  data            <- mon_inits$data
  extension_inits <- mon_inits$extension_inits


  # Store MCMC compilation and run results
  samples <- NULL
  model <- NULL

  nimble_waic <- NULL
  if (engine == "nimble") {
    # --- NIMBLE EXECUTION PIPELINE ---
    nimble_res <- run_nimble_pipeline(
      model_string    = model_string,
      data            = data,
      family          = family,
      extension_inits = extension_inits,
      monitor         = monitor,
      n.chains        = n.chains,
      n.iter          = n.iter,
      n.burnin        = n.burnin,
      n.thin          = n.thin,
      WAIC            = WAIC,
      nimble_samplers = nimble_samplers,
      parallel        = parallel,
      n.cores         = n.cores,
      cl              = cl,
      quiet           = quiet
    )
    samples               <- nimble_res$samples
    model                 <- nimble_res$model
    monitor               <- nimble_res$monitor
    nimble_waic           <- nimble_res$nimble_waic
    saved_nimble_compiled <- nimble_res$saved_nimble_compiled
    saved_nimble_cmodel   <- nimble_res$saved_nimble_cmodel
    saved_nimble_samplers <- nimble_res$saved_nimble_samplers
  } else {
    # --- JAGS EXECUTION PIPELINE (Default) ---
    jags_res <- run_jags_pipeline(
      model_file = model_file,
      model_string = model_string,
      data = data,
      extension_inits = extension_inits,
      monitor = monitor,
      n.chains = n.chains,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.thin = n.thin,
      n.adapt = n.adapt,
      quiet = quiet,
      verbose = verbose,
      parallel = parallel,
      n.cores = n.cores,
      cl = cl,
      DIC = DIC,
      WAIC = WAIC
    )
    samples <- jags_res$samples
    model   <- jags_res$model
    DIC     <- jags_res$DIC
    WAIC    <- jags_res$WAIC
  }

  return(assemble_because_result(
    model                 = model,
    model_code            = model_string,
    model_file            = model_file,
    samples               = samples,
    data                  = data,
    original_data         = original_data,
    equations             = equations,
    random                = random,
    random_terms          = random_terms,
    structure             = structure,
    structures            = structures,
    latent                = latent,
    distribution          = distribution,
    family                = family,
    variability           = variability,
    all_poly_terms        = all_poly_terms,
    dsep                  = dsep,
    dsep_tests            = dsep_tests,
    dsep_results          = dsep_results,
    parameter_map         = parameter_map,
    induced_cors          = induced_cors,
    scale_info            = scale_info,
    stack_res             = if (exists("stack_res")) stack_res else NULL,
    engine                = engine,
    saved_nimble_compiled = if (exists("saved_nimble_compiled")) saved_nimble_compiled else NULL,
    saved_nimble_cmodel   = if (exists("saved_nimble_cmodel")) saved_nimble_cmodel else NULL,
    saved_nimble_samplers = if (exists("saved_nimble_samplers")) saved_nimble_samplers else NULL,
    nimble_waic           = if (exists("nimble_waic")) nimble_waic else NULL,
    parallel              = parallel,
    n.cores               = n.cores,
    n.chains              = n.chains,
    n.iter                = n.iter,
    n.burnin              = n.burnin,
    n.thin                = n.thin,
    n.adapt               = n.adapt,
    DIC                   = DIC,
    WAIC                  = WAIC,
    ic_recompile          = ic_recompile,
    extension_inits       = extension_inits,
    quiet                 = quiet,
    id_col                = id_col,
    original_call         = original_call,
    hierarchical_info     = hierarchical_info,
    monitor               = monitor
  ))
}

