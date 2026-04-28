#' Harden NIMBLE Sampler Configuration
#'
#' Applies robust sampler assignments to a NIMBLE MCMC configuration object.
#' Exported as a standalone function so that parallel worker nodes — which
#' load the package fresh — always use the current installed version of this
#' logic, regardless of which version of `because()` was originally called.
#'
#' @param mcmc_conf A NIMBLE MCMC configuration object (from `configureMCMC()`).
#' @param family Named character vector of response families (same as `because()`).
#' @param nimble_samplers Optional named list of user-specified samplers.
#' @param quiet Logical. If TRUE, suppress status messages.
#' @return The modified `mcmc_conf` object (invisibly).
#' @keywords internal
#' @export
nimble_harden_samplers <- function(mcmc_conf, family = NULL, nimble_samplers = NULL, quiet = TRUE) {

  # Helper: parse VAR name from an err_raw_VAR_STRUCTURE[...] or u_std_VAR_STRUCTURE[...] node
  # Strategy: strip prefix and index, then the last _-delimited token is the structure name.
  parse_re_info <- function(node) {
    base   <- sub("\\[.*\\]", "", node)                  # strip [1:N]
    base   <- sub("^(err_raw_|u_std_|sigma_|tau_|beta_|alpha_)", "", base) # strip prefixes
    parts  <- strsplit(base, "_")[[1]]
    if (length(parts) < 2) return(list(var = base, structure = ""))
    
    # For betas: beta_RESPONSE_PREDICTOR
    if (grepl("^beta_", node)) {
       return(list(var = parts[1], structure = "beta"))
    }

    # For REs and scales: name_STRUCTURE
    last   <- parts[length(parts)]
    known_structs <- c("phylo", "spatial", "survey", "site", "obs", "res")
    if (last %in% known_structs) {
       return(list(
         var       = paste(parts[-length(parts)], collapse = "_"), 
         structure = last                                          
       ))
    }
    list(var = base, structure = "")
  }

  sampler_targets <- sapply(mcmc_conf$getSamplers(), function(x) x$target)

  # ── 1. Hybrid Grouping Logic: Core Equation Blocks ────────────────────────
  # The goal is a SMALL block (5-10 nodes) containing: (alpha, all betas, all scales)
  trait_groups <- list()

  # Identify possible intercepts
  alpha_nodes <- unique(grep("^alpha_.*", sampler_targets, value = TRUE))
  for (a_node in alpha_nodes) {
    trait_name <- sub("^alpha_", "", a_node)
    trait_groups[[trait_name]] <- list(targets = a_node)
  }

  # Identify all Fixed Effects (Slopes) and Variance nodes
  hyper_nodes <- unique(grep("^(beta_|sigma_|tau_).*", sampler_targets, value = TRUE))
  
  for (node in hyper_nodes) {
    # Skip if already handled by posterior_predictive logic
    current_types <- sapply(mcmc_conf$getSamplers(node), function(s) s$name)
    if (any(grepl("posterior_predictive", current_types, ignore.case = TRUE))) next

    info <- parse_re_info(node)
    t_name <- info$var
    
    # Heuristic for beta_RESPONSE_PREDICTOR: find the matching trait name
    if (grepl("^beta_", node)) {
       # Find which trait_group it belongs to
       matches <- names(trait_groups)[vapply(names(trait_groups), function(tn) grepl(paste0("^", tn, "_"), sub("^beta_", "", node)), logical(1))]
       if (length(matches) > 0) t_name <- matches[which.max(nchar(matches))]
    }

    # Clean up standard trait name prefixing
    if (t_name == "" || is.na(t_name)) next
    if (grepl("^phylo_", t_name)) t_name <- sub("^phylo_", "", t_name)

    # Add to group if a matching intercept (alpha) exists
    if (!is.null(trait_groups[[t_name]])) {
       trait_groups[[t_name]]$targets <- unique(c(trait_groups[[t_name]]$targets, node))
    }
  }

  # ── 2. Apply Joint Core Blocks (AF_slice vs RW_block) ─────────────────────────
  processed_nodes <- character(0)

  for (t_name in names(trait_groups)) {
    group_targets <- trait_groups[[t_name]]$targets
    
    # We block if there's an intercept AND at least one scale/slope component
    if (length(group_targets) > 1) {
      for (target in group_targets) {
        mcmc_conf$removeSamplers(target)
      }

      # Root Node Detection: Traits with no slopes (betas) get AF_slice
      has_slopes   <- any(grepl("^beta_", group_targets))
      sampler_type <- if (has_slopes) "RW_block" else "AF_slice"

      mcmc_conf$addSampler(
        target = group_targets,
        type   = sampler_type
      )
      
      processed_nodes <- c(processed_nodes, group_targets)

      if (!quiet) {
        message(sprintf(
          "NIMBLE: %s Core Block for '%s' (%d nodes: %s, etc.)",
          sampler_type, t_name, length(group_targets), group_targets[1]
        ))
      }
    }
  }

  # ── 3. Handle Remaining Nodes (ESS and Defaults) ──────────────────────────
  remaining_targets <- setdiff(sampler_targets, processed_nodes)
  
  for (target in remaining_targets) {
    current_types <- sapply(mcmc_conf$getSamplers(target), function(s) s$name)
    if (any(grepl("posterior_predictive", current_types, ignore.case = TRUE))) next
    
    # [NEW] Phylogenetic ESS: Force Elliptical Slice Sampler for multivariate phylo nodes
    if (grepl("^(err_raw_|u_std_).*(phylo|BM|OU|Pagel).*", target)) {
        mcmc_conf$removeSamplers(target)
        mcmc_conf$addSampler(target = target, type = "ess")
        if (!quiet) message(sprintf("NIMBLE: Forced 'ess' (Elliptical Slice) for phylogenetic vector '%s'.", target))
        next
    }

    # Other REs (Standalone / Default)
    if (grepl("^(err_raw_|u_std_).*", target)) {
        if (length(current_types) > 1 || (!any(grepl("RW_block|conjugate|ess", current_types, ignore.case = TRUE)))) {
            mcmc_conf$removeSamplers(target)
            mcmc_conf$addSampler(target = target, type = "RW_block")
        }
    }
    
    # Standalone Scales (Standalone Defaults)
    if (grepl("^(sigma_|tau_|lambda_|r_|psi_|sigma_total_).*", target)) {
      mcmc_conf$removeSamplers(target)
      mcmc_conf$addSampler(target = target, type = "slice")
    }
  }

  if (!quiet) {
      message(sprintf("NIMBLE: Sampler hardening complete. ESS + AF_slice + Hybrid strategy applied."))
  }

  # ── 7. User-specified overrides (always applied last) ─────────────────────
  if (!is.null(nimble_samplers)) {
    for (node in names(nimble_samplers)) {
      mcmc_conf$removeSamplers(node)
      mcmc_conf$addSampler(target = node, type = nimble_samplers[[node]])
    }
  }

  invisible(mcmc_conf)
}


#' @title Run a Bayesian Structural Equation Model (Because)
#'
#' @description
#' Fits a Bayesian Structural Equation Model (SEM) using JAGS or NIMBLE.
#' Supports multi-level (hierarchical) data, custom covariance structures
#' (phylogenetic, spatial, etc.), missing data imputation, and d-separation
#' global fit testing.
#'
#' @param equations A list of R formulas describing the structural model. Each formula represents 
#'    a causal path (e.g., `Y ~ X1 + X2`). Categorical predictors are automatically handled.
#' @param data A data frame containing all variables. For **hierarchical models**, 
#'   this must be a named list of data frames (one per level).
#' @param family A named character vector specifying the distribution for each response 
#'   (e.g., `c(Y = "poisson", M = "gaussian")`). Supports "gaussian" (default), "poisson", 
#'   "binomial", "gamma", "lognormal", "bernoulli", "ordinal", and "occupancy".
#' @param dsep Logical; if `TRUE`, performs d-separation tests to evaluate the global 
#'   fit of the DAG against the data. Highly recommended for causal validation.
#' @param random Optional formula for **global random intercepts** (e.g., `~(1|Site)`). 
#'   Deprecated in favor of inline specification: `Y ~ X + (1|Site)`.
#' @param latent Optional character vector of latent (unmeasured) variables.
#' @param id_col Character string for the column identifying units (e.g. species names).
#' @param structure Optional covariance structure (e.g., a phylogenetic tree or spatial matrix) 
#'   for modeling correlated residuals.
#' @param hierarchy (Hierarchical) A string defining the nesting structure (e.g., `"site > individual"`). 
#'   Semicolons separate parallel branches (e.g., `"site > obs; species > obs"`).
#' @param levels (Hierarchical) Optional named list mapping variables to their home levels. 
#'   If omitted, `because` will attempt to auto-detect levels based on column availability.
#' @param link_vars (Hierarchical) A named character vector specifying the columns used to 
#'   link data frames across levels (e.g., `c(site = "SiteID")`).
#' @param engine Bayesian backend: `"jags"` (default) or `"nimble"`. 
#' @param n.iter Total MCMC iterations per chain (default = 12500).
#' @param n.burnin Number of burn-in iterations (default = 20% of `n.iter`).
#' @param n.chains Number of independent MCMC chains (default = 3).
#' @param parallel Logical; if `TRUE`, runs MCMC chains in parallel.
#' @param n.cores Number of CPU cores for parallel execution.
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
#' # 2. Hierarchical Model (Auto-detection)
#' # Suppose we have year-level data and individual-level data
#' year_df <- data.frame(Year = 1:5, Temp = rnorm(5))
#' ind_df  <- data.frame(Year = rep(1:5, each=10), Mass = rnorm(50))
#' data_list <- list(yr = year_df, ind = ind_df)
#' 
#' # Mass depends on Temp (cross-level link via Year)
#' fit_h <- because(
#'   equations = list(Mass ~ Temp),
#'   data      = data_list,
#'   hierarchy = "yr > ind",
#'   link_vars = c(yr = "Year")
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
#' @importFrom utils capture.output head
#' @importFrom coda gelman.diag effectiveSize
because <- function(
  equations,
  data,
  id_col = NULL,
  structure = NULL,
  engine = "jags",
  monitor = "interpretable",
  nimble_samplers = NULL,
  n.chains = 3,
  n.iter = 12500,
  n.burnin = floor(n.iter / 5),
  n.thin = 10,
  DIC = TRUE,
  WAIC = FALSE,
  n.adapt = floor(n.iter / 5),
  quiet = FALSE,
  verbose = FALSE,
  dsep = FALSE,
  variability = NULL,
  family = NULL,
  distribution = NULL,
  latent = NULL,
  latent_method = "correlations",
  standardize_latent = TRUE,
  fix_latent = "loading",
  parallel = FALSE,
  n.cores = parallel::detectCores() - 1,
  cl = NULL,
  ic_recompile = FALSE,
  random = NULL,
  levels = NULL,
  hierarchy = NULL,
  link_vars = NULL,
  fix_residual_variance = NULL,
  priors = NULL,
  reuse_models = NULL,
  expand_ordered = FALSE,
  structure_multi = NULL,
  structure_levels = NULL,
  ...
) {
  args <- list(...)
  
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

  engine <- match.arg(tolower(engine), c("jags", "nimble"))

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

  # --- Clean Data: Ensure matrix (e.g. from scale()) is a data.frame ---
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # --- Clean Data: Coerce 1D matrices (e.g., from scale()) to vectors ---
  if (is.data.frame(data)) {
    for (nm in names(data)) {
      if (is.matrix(data[[nm]]) && ncol(data[[nm]]) == 1) {
        data[[nm]] <- as.numeric(data[[nm]])
      }
    }
  }

  # --- Backward Compatibility: distribution -> family ---
  if (!is.null(distribution)) {
    warning(
      "Argument 'distribution' is deprecated and will be removed in future versions. Please use 'family' instead."
    )
    if (is.null(family)) {
      family <- distribution
    }
  }

  # --- Family Object Normalization (Custom Families Support) ---
  # Store original family objects for S3 dispatch, normalize to character vector
  family_objects <- list()
  if (!is.null(family)) {
    if (inherits(family, "because_family")) {
      # Single family object for all responses
      family_objects[["_default"]] <- family
      family <- setNames(family$family, "_default")
    } else if (is.list(family) && !is.null(names(family))) {
      # Named list: check each element for family objects
      for (nm in names(family)) {
        if (inherits(family[[nm]], "because_family")) {
          family_objects[[nm]] <- family[[nm]]
          family[[nm]] <- family[[nm]]$family
        }
      }
      family <- unlist(family)
    }
    # Now family is a named character vector, family_objects stores originals
  }

  # S3 Class assignment for dispatch on extensions
  family_obj <- family
  if (!is.null(family) && ("occupancy" %in% family || "cmr" %in% family)) {
    class(family_obj) <- c("because_family_occupancy", class(family_obj))
  }
  structure_obj <- structure
  if (!is.null(structure)) {
    class(structure_obj) <- c("because_structure", class(structure_obj))
  }

  # Extension Hook: Normalize specialized aliases (e.g. psi_Species -> Species)
  equations <- normalize_equations_hook(family_obj, equations, data = data)

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
        data <- data[available_vars]
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


      # Use auto-detected hierarchy/link_vars if not explicitly provided
      if (is.null(hierarchy)) {
        hierarchy <- auto_result$hierarchy
      }
      if (is.null(link_vars)) {
        link_vars <- auto_result$link_vars
      }
    }

    # Now validate (with either provided or auto-detected values)
    if (!is.null(levels) && length(levels) > 0) {
      is_hierarchical <- TRUE
      validate_hierarchical_data(
        data,
        levels,
        hierarchy,
        link_vars,
        latent_vars = latent,
        equations = equations
      )
      # Try to infer hierarchy from random effects if still not set
      if (is.null(hierarchy)) {
        hierarchy <- parse_hierarchy_from_random(random, data)
        if (is.null(hierarchy)) {
          stop(
            "Hierarchical data detected but 'hierarchy' not specified. ",
            "Provide either:\n",
            "  1. 'hierarchy' argument (e.g., \"site_year > individual\"), or\n",
            "  2. Nested random effects (e.g., ~(1|site/individual))"
          )
        }
      }

      # Store hierarchical info for later use
      hierarchical_info <- list(
        data = data,
        levels = levels,
        hierarchy = hierarchy,
        link_vars = link_vars
      )

      # Internal: Inject structural metadata if provided in sub-calls (e.g. from d-sep tests)
      if (!is.null(structure_multi)) {
        hierarchical_info$structure_multi <- structure_multi
      }
      if (!is.null(structure_levels)) {
        hierarchical_info$structure_levels <- structure_levels
      }

      if (!quiet) {
        message("Hierarchical data structure detected: ", hierarchy)
      }
    }
  } else {
    # Single-level data - check if hierarchical metadata was provided anyway (e.g. for d-sep tests)
    if (!is.null(levels) && !is.null(hierarchy)) {
      is_hierarchical <- FALSE # Data is flat, so NOT hierarchical for prep purposes
      hierarchical_info <- list(
        data = data,
        levels = levels,
        hierarchy = hierarchy,
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
      "See vignette('07_multilevel_models') for details.",
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
        "Prepared hierarchical data for JAGS: ",
        paste(
          names(prep_res$n_vec),
          unlist(prep_res$n_vec),
          sep = "=",
          collapse = ", "
        )
      )
    }

    # original_data remains the list of dataframes for safety
    original_data <- hierarchical_info$data
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

  # --- Random Effects Data Prep (Post-Assembly) ---
  # Create structures for JAGS using the assembled data
  if (length(random_terms) > 0) {
    rand_structs <- create_group_structures(data, random_terms)

    if (!quiet) {}
    rand_structs <- create_group_structures(data, random_terms)
    random_structures <- rand_structs$structures
    random_data_updates <- rand_structs$data_updates
  }

  if (is.data.frame(data) || (is.list(data) && !is.data.frame(data))) {
    # Check for long format data requiring matrix conversion
    # If any variability specified as 'reps', we attempt to auto-format using because_format_data
    has_reps <- any(grepl("reps", variability))

    if (has_reps) {
      if (is.null(structure)) {
        # Cannot auto-format without structure to determine species order
        warning(
          "Variability 'reps' specified but no structure provided. Automatic formatting requires a structure to order species rows. Assuming data is already aggregated or user handles index mapping."
        )
        return(data)
      } else if (!is.null(id_col) && id_col %in% names(data)) {
        # Extension Hook: Extract appropriate tree from structure object
        use_tree <- get_tree_hook(structure)

        formatted_list <- because_format_data(
          data,
          species_col = id_col,
          tree = use_tree
        )

        # Now rename the variables that are 'reps' to include '_obs' suffix
        # And keep others as is
        reps_vars <- names(variability)[variability == "reps"]

        final_data_list <- list()
        for (nm in names(formatted_list)) {
          if (nm %in% reps_vars) {
            # This is a matrix of replicates, rename to _obs
            final_data_list[[paste0(nm, "_obs")]] <- formatted_list[[nm]]
          } else {
            # This is a regular variable (vector or matrix depending on because_format_data logic)
            final_data_list[[nm]] <- formatted_list[[nm]]
          }
        }

        # Update data to be this list
        data <- final_data_list

        if (!quiet) {
          message(
            "  Formatted ",
            length(names(formatted_list)),
            " variables as replicate matrices."
          )
        }

        # Check for missing variables that were expected to be formatted
        # Filter variability/family to variables actually in current equations
        # to avoid 'missing reps' errors for irrelevant variables in sub-models (e.g. d-sep tests)
        all_vars <- unique(c(
          names(equations),
          unlist(lapply(equations, all.vars))
        ))
        # Include versions with p_ removed for variability/dist matching
        clean_vars <- sub("^p_", "", all_vars)
        relevant_vars <- unique(c(all_vars, clean_vars))

        missing_reps <- setdiff(reps_vars, names(formatted_list))
        # Only error if the missing variable is actually relevant to our equations
        missing_reps <- intersect(missing_reps, relevant_vars)

        if (length(missing_reps) > 0) {
          stop(paste(
            "The following variables were identified for 'reps' processing (from equations) but were NOT found in the data:",
            paste(missing_reps, collapse = ", "),
            "\nPlease check your column names."
          ))
        }
      } else {
        warning(
          "Variability 'reps' specified but 'id_col' missing or not in data. Cannot auto-format long data."
        )
      }
    }
  }

  if (
    (is.data.frame(data) || (is.list(data) && !is.data.frame(data))) &&
      !is_hierarchical
  ) {
    # Extract all variable names from fixed equations
    eq_vars <- unique(unlist(lapply(equations, all.vars)))

    # Add variables from random terms (grouping factors)
    if (length(random_terms) > 0) {
      random_vars <- unique(vapply(
        random_terms,
        function(x) x$group,
        character(1)
      ))
      eq_vars <- unique(c(eq_vars, random_vars))
    }

    # Check which variables are in the data frame
    # Check which variables are in the data frame
    # We must include "N" if it exists, as it's needed for optimized loops
    cols_to_keep <- unique(c(eq_vars, "N"))
    if (!is.null(variability) && !is.character(variability)) {
      cols_to_keep <- c(cols_to_keep, names(variability))
    }

    available_vars <- intersect(
      cols_to_keep,
      (if (is.null(names(data))) character(0) else names(data))
    )
    missing_vars <- setdiff(
      eq_vars,
      (if (is.null(names(data))) character(0) else names(data))
    )

    # Some "missing" vars might be latent - that's OK
    if (!is.null(latent)) {
      missing_vars <- setdiff(missing_vars, latent)
    }

    if (length(missing_vars) > 0 && length(available_vars) == 0) {
      stop(
        "None of the variables in equations found in data frame. ",
        "Missing: ",
        paste(missing_vars, collapse = ", ")
      )
    }

    if (length(missing_vars) > 0 && !quiet) {
      message(
        "Note: Variables not in data (may be latent/derived): ",
        paste(missing_vars, collapse = ", ")
      )
    }

    # Handle id_col for matching to tree/structure
    row_ids <- NULL
    if (is.data.frame(data)) {
      if (!is.null(id_col)) {
        if (!id_col %in% names(data)) {
          stop("id_col '", id_col, "' not found in data frame columns.")
        }
        row_ids <- data[[id_col]]
        # Remove id_col from variables to include (it's metadata, not a model variable)
        available_vars <- setdiff(available_vars, id_col)
      } else {
        # Try to use row names if they're meaningful (not just 1, 2, 3...)
        rn <- rownames(data)
        if (!is.null(rn) && !all(rn == as.character(seq_len(nrow(data))))) {
          row_ids <- rn
        }
      }
    }
    # Also check for variability-related columns (X_se, X_obs patterns)
    se_cols <- grep("_se$", names(data), value = TRUE)
    obs_cols <- grep("_obs$", names(data), value = TRUE)

    # Add generated dummy variables for categorical predictors
    dummy_vars <- character(0)
    if (!is.null(attr(data, "categorical_vars"))) {
      cat_vars <- attr(data, "categorical_vars")

      # Extract RHS variables (predictors) from fixed equations to filter dummies
      rhs_vars <- unique(unlist(lapply(equations, function(eq) {
        if (length(eq) == 3) {
          all.vars(eq[[3]])
        } else {
          character(0)
        }
      })))

      # Only generate dummies for variables used as predictors
      cat_vars <- cat_vars[names(cat_vars) %in% rhs_vars]

      dummy_vars <- unlist(lapply(cat_vars, function(x) x$dummies))
    }

    extra_cols <- c(se_cols, obs_cols, dummy_vars)

    # Convert to list format

    data_list <- list()
    for (var in c(available_vars, extra_cols)) {
      if (var %in% names(original_data)) {
        data_list[[var]] <- original_data[[var]]
      } else if (var %in% names(data)) {
        data_list[[var]] <- data[[var]]
      }
    }

    # Set names on vectors if we have row_ids
    if (!is.null(row_ids)) {
      for (var in names(data_list)) {
        if (
          is.vector(data_list[[var]]) &&
            length(data_list[[var]]) == length(row_ids)
        ) {
          names(data_list[[var]]) <- row_ids
        }
      }
    }

    # Preserve any existing attributes
    data_attrs <- attributes(original_data)

    # Combine data_list with random effect data
    data <- data_list
    if (length(random_data_updates) > 0) {
      for (nm in names(random_data_updates)) {
        if (!is.null(nm) && nm != "") {
          data[[nm]] <- random_data_updates[[nm]]
        }
      }
    }

    # Remove raw grouping variables from data passed to JAGS to avoid "Unused variable" warnings
    if (length(random_terms) > 0) {
      # Use setdiff to avoid error if variable not present (though they should be)
      vars_to_remove <- intersect(names(data), random_vars)
      if (length(vars_to_remove) > 0) {
        data[vars_to_remove] <- NULL
      }
    }

    if (!quiet) {
      message(
        "Converted data.frame to list with ",
        length(data),
        " variables: ",
        paste(names(data), collapse = ", ")
      )
    }
  } else {
    # If optimized hierarchical path was taken, we skipped the big block above.
    # We MUST merge random data updates (Prec matrices) now.
    if (is_hierarchical) {
      if (length(random_data_updates) > 0) {
        for (nm in names(random_data_updates)) {
          if (!is.null(nm) && nm != "") {
            data[[nm]] <- random_data_updates[[nm]]
          }
        }
      }
      # Restore attributes (categorical_vars) needed for model generation
      if (
        !is.null(hierarchical_info) &&
          !is.null(attr(hierarchical_info$data, "categorical_vars"))
      ) {
        attr(data, "categorical_vars") <- attr(
          hierarchical_info$data,
          "categorical_vars"
        )
      }
      data_attrs <- list()
    } else {
      # Standard Fallback
      # Ensure data is a list (crucial for adding matrices like VCV)
      # Preserve attributes (like categorical_vars) which are lost during as.list()
      data_attrs <- attributes(data)
      data <- as.list(data)
    }
  }

  # Compute polynomial values
  # We add them to original_data (for d-sep/residuals) but NOT to 'data' passed to JAGS
  # because JAGS creates deterministic nodes for them (var_pow2 <- var^2)
  if (!is.null(all_poly_terms)) {
    for (poly_term in all_poly_terms) {
      base_var <- poly_term$base_var
      power <- poly_term$power
      internal_name <- poly_term$internal_name

      # Check if base variable exists in data
      if (base_var %in% names(data)) {
        # Compute polynomial: x^2, x^3, etc.
        poly_vals <- data[[base_var]]^power

        # Add to original_data if possible
        if (is.list(original_data) || is.data.frame(original_data)) {
          original_data[[internal_name]] <- poly_vals
        }

        # If hierarchical, we MUST also add it to the source dataframes in hierarchical_info
        # Otherwise d-separation tests (which re-fetch data) won't find the new variable
        if (is_hierarchical && !is.null(hierarchical_info)) {
          # Find which level/dataframe holds the base variable
          for (lvl_name in names(hierarchical_info$data)) {
            if (base_var %in% names(hierarchical_info$data[[lvl_name]])) {
              # Add computed column to this level's dataframe
              source_df <- hierarchical_info$data[[lvl_name]]
              hierarchical_info$data[[lvl_name]][[internal_name]] <- source_df[[
                base_var
              ]]^power
              break
            }
          }
        }
      }
    }
  }

  # Restore categorical_vars if present
  if ("categorical_vars" %in% names(data_attrs)) {
    attr(data, "categorical_vars") <- data_attrs$categorical_vars
  }

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
      if (is.list(s) && length(s) > 1 && !is.matrix(s) && !inherits(s, "phylo")) {
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
    structure_levels <- list()
    structure_multi <- list()
    # Use S3 Generic for Processing
    for (s_name in structure_names) {
      structure_obj <- structures[[s_name]]
      prep_res <- prepare_structure_data(structure_obj, data = data, optimize = TRUE, quiet = quiet)

      if (!is.null(prep_res$data_list)) {
        for (d_name in names(prep_res$data_list)) {
          custom_s_name <- get_structure_name_hook(structure_obj)
          prefixed_name <- if (d_name %in% c("Prec", "VCV", "multiVCV", custom_s_name)) paste0(d_name, "_", s_name) else d_name
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
      hierarchical_info$structure_levels <- structure_levels
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
    if (var %in% names(data)) {
      var_data <- data[[var]]
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

  # Handle d-sep logic
  dsep_tests <- NULL
  induced_cors <- NULL

  if (dsep) {
    # Extension Hook: Remove specialized variables from potential latents

    # Force WAIC and DIC off for d-separation testing (not needed for conditional independence tests)
    if (WAIC || DIC) {
      if (!quiet) {
        message(
          "Note: WAIC and DIC are not computed for d-separation tests (not needed for conditional independence testing)."
        )
      }
      WAIC <- FALSE
      DIC <- FALSE
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
          message(
            "Auto-detected latent variable(s): ",
            paste(latent, collapse = ", "),
            "\n(Variables in equations but not in data will be treated as latent.)\n",
            "Generating m-separation tests for MAG..."
          )
        }
      }
    }

    if (!quiet) {
      if (!is.null(latent)) {
        message(
          "Generating m-separation tests (MAG with latent variables)..."
        )
      } else {
        message("Generating d-separation tests...")
      }
    }
    # Expanded random terms for d-sep (include all variables as potential responses)
    # This ensures that root nodes in the DAG get random effects in d-sep tests
    random_terms_for_dsep <- random_terms
    if (!is.null(random)) {
      # Use all variables currently in the model equations
      all_vars_in_model <- unique(unlist(lapply(equations, all.vars)))
      extra_random <- parse_global_random(
        random,
        equations,
        all_vars = all_vars_in_model
      )

      # Combine and deduplicate
      combined_rand <- c(random_terms, extra_random)
      rand_keys <- sapply(combined_rand, function(x) {
        paste(x$response, x$group, sep = "|")
      })
      random_terms_for_dsep <- combined_rand[!duplicated(rand_keys)]
    }

    dsep_result <- because_dsep(
      equations,
      latent = latent,
      random_terms = random_terms_for_dsep,
      categorical_vars = if (!is.null(attr(data, "categorical_vars"))) {
        attr(data, "categorical_vars")
      } else {
        NULL
      },
      family = family,
      poly_terms = all_poly_terms,
      quiet = quiet,
      hierarchical_info = hierarchical_info
    )

    # Extract tests and correlations
    if (!is.null(latent)) {
      dsep_tests <- dsep_result$tests
      induced_cors <- dsep_result$correlations

      if (!quiet && length(induced_cors) > 0) {
        message(
          "Found ",
          length(induced_cors),
          " induced correlation(s) from latent variable(s)"
        )
      }
    } else {
      dsep_tests <- dsep_result
    }

    # Extension Hook: Translate tests where response is specialized (e.g. psi_)
    dsep_tests <- lapply(dsep_tests, function(eq) {
      attr_val <- attr(eq, "test_var")
      new_eq <- dsep_test_translation_hook(family_obj, eq)
      if (!is.null(attr_val)) attr(new_eq, "test_var") <- attr_val
      new_eq
    })

    # (Original translation block removed)

    # Deduplicate tests in case translation created duplicates
    dsep_test_strs <- sapply(dsep_tests, function(eq) {
      # Use trimws and a single space to normalize for string comparison
      str <- paste(deparse(eq), collapse = " ")
      gsub("\\s+", " ", trimws(str))
    })

    dsep_tests <- dsep_tests[!duplicated(dsep_test_strs)]

    if (length(dsep_tests) == 0) {
      stop(
        "No d-separation tests implied by the model (model is saturated). Stopping run."
      )
    }

    # --- A priori cross-hierarchy filter ---
    # Tests where the response and the focal predictor belong to orthogonal
    # hierarchical branches (e.g., species-level vs. site-level) are trivially
    # satisfied by construction: there is no causal path between the two branches
    # except through the observation level, and no cross-level index exists in
    # the data to estimate such a regression. Skip them with a clear message
    # rather than letting JAGS fail with "Unknown variable species_idx_site".
    cross_hier_flags <- sapply(dsep_tests, function(eq) {
      is_cross_hierarchy_test(eq, hierarchical_info)
    })
    n_cross <- sum(cross_hier_flags)
    if (n_cross > 0 && !quiet) {
      skipped_labels <- sapply(dsep_tests[cross_hier_flags], function(eq) {
        resp     <- as.character(eq)[2]
        test_var <- attr(eq, "test_var")
        if (!is.null(test_var)) {
          sprintf("  %s _||_ %s (orthogonal hierarchy branches - trivially satisfied)",
                  resp, test_var)
        } else {
          paste(deparse(eq), collapse = " ")
        }
      })
      message(sprintf(
        "\nSkipping %d cross-hierarchy d-sep test(s) (trivially satisfied by design):",
        n_cross
      ))
      for (lbl in skipped_labels) message(lbl)
    }
    dsep_tests_skipped <- dsep_tests[cross_hier_flags]
    dsep_tests <- dsep_tests[!cross_hier_flags]

    if (length(dsep_tests) == 0) {
      warning("All d-separation tests were cross-hierarchy (trivially satisfied). ",
              "No JAGS models will be run.")
    }

    # Run tests sequentially to avoid cyclic dependencies in JAGS
    # Decide whether to run tests in parallel
    use_parallel <- parallel && n.cores > 1 && length(dsep_tests) > 1

    # INCREMENTAL D-SEP: Check for reusable tests
    reused_results <- list()
    tests_to_run_indices <- seq_along(dsep_tests)

    if (!is.null(reuse_models)) {
      incremental_check <- find_reusable_tests(
        dsep_tests,
        equations,
        reuse_models,
        data,
        family = family,
        quiet = quiet
      )
      reused_results <- incremental_check$found # List with results at matching indices, NULL otherwise
      tests_to_run_indices <- incremental_check$missing_indices

      # Override parallel decision if few tests remain
      if (length(tests_to_run_indices) < 2) {
        use_parallel <- FALSE
      }
    }

    if (!quiet) {
      if (use_parallel) {
        message(sprintf(
          "Running %d d-separation tests in parallel on %d cores...",
          length(tests_to_run_indices), # Corrected message to show actual runs
          n.cores
        ))
      } else {
        message(sprintf(
          "Running %d d-separation tests sequentially...",
          length(tests_to_run_indices)
        ))
      }
    }

    combined_samples <- NULL
    combined_map <- NULL

    # Define function to run a single d-sep test
    # run_single_dsep_test_v2 removed from here (now top-level)

    # Run tests (parallel or sequential)
    # Only loop over tests_to_run_indices
    # We still need to populate a full results list of length(dsep_tests)
    # The 'reused_results' list has NULLs for missing tests and values for found ones.

    new_results_list <- vector("list", length(dsep_tests))

    # Fill in reused results first
    if (length(reused_results) > 0) {
      for (i in seq_along(reused_results)) {
        if (!is.null(reused_results[[i]])) {
          new_results_list[[i]] <- reused_results[[i]]
        }
      }
    }

    if (length(tests_to_run_indices) > 0) {
      if (use_parallel) {
        # Setup cluster if not provided
        if (is.null(cl)) {
          cl <- parallel::makeCluster(n.cores)
          on.exit(parallel::stopCluster(cl), add = TRUE)
        }

        # Identify because extensions loaded in the master session
        parent_exts <- grep("^because\\.", loadedNamespaces(), value = TRUE)

        # Ensure workers have the library loaded
        parallel::clusterEvalQ(cl, {
          if (requireNamespace("because", quietly = TRUE)) {
            library(because)
          }
        })

        # Load detected extensions on workers
        if (length(parent_exts) > 0) {
            # Pass the list of extensions to workers
            parallel::clusterExport(cl, "parent_exts", envir = environment())
            parallel::clusterEvalQ(cl, {
                for (ext in parent_exts) {
                    if (requireNamespace(ext, quietly = TRUE)) {
                        library(ext, character.only = TRUE)
                    }
                }
            })
        }

        # Export necessary objects to cluster
        parallel::clusterExport(
          cl,
          c(
            "original_data",
            "structure",
            "monitor",
            "n.chains",
            "n.iter",
            "n.burnin",
            "n.thin",
            "n.adapt",
            "variability",
            "family",
            "latent",
            "latent_method",
            "ic_recompile",
            "equations",
            "random_terms",
            "hierarchical_info",
            "levels",
            "hierarchy",
            "link_vars",
            "fix_residual_variance",
            "run_single_dsep_test_v2",
            "dsep_tests",
            "extract_random_effects",
            "engine",
            "nimble_samplers",
            "because",
            "quiet",
            "random",
            "get_data_for_variables",
            "infer_variable_level",
            "get_level_depth",
            "sanitize_term_name",
            "dsep_tree_hook",
            "dsep_equations_hook"
          ),
          envir = environment()
        )

        # Run tests in parallel
        # Only run the necessary tests
        # We use pbapply for a live progress bar if available, otherwise fallback to parLapply
        # Note: pblapply and parLapply have different argument orders for 'cl'
        if (requireNamespace("pbapply", quietly = TRUE)) {
          par_results <- pbapply::pblapply(
            X = tests_to_run_indices,
            cl = cl,
            FUN = function(i) {
            test_eq <- dsep_tests[[i]]

            # Use "interpretable" monitoring mode.
            current_monitor <- "interpretable"

            tryCatch({
              run_single_dsep_test_v2(
                i,
                test_eq,
                current_monitor,
                engine = engine,
                nimble_samplers = nimble_samplers,
                quiet = quiet,
                original_data = original_data,
                hierarchical_info = hierarchical_info,
                random_terms = random_terms,
                equations = equations,
                family = family,
                structure = structure,
                levels = levels,
                hierarchy = hierarchy,
                link_vars = link_vars,
                fix_residual_variance = fix_residual_variance,
                latent = latent,
                latent_method = latent_method,
                n.chains = n.chains,
                n.iter = n.iter,
                n.burnin = n.burnin,
                n.thin = n.thin,
                n.adapt = n.adapt,
                ic_recompile = ic_recompile,
                random = random,
                id_col = id_col,
                variability = variability
              )
            }, error = function(e) {
              # Return error info to parent session for reporting
              list(
                error = conditionMessage(e),
                test_idx = i,
                equation = paste(deparse(test_eq), collapse = " ")
              )
            })
          }
        )
        } else {
          par_results <- parallel::parLapply(
            cl = cl,
            X = tests_to_run_indices,
            fun = function(i) {
              test_eq <- dsep_tests[[i]]

              # Use "interpretable" monitoring mode.
              current_monitor <- "interpretable"

              tryCatch({
                run_single_dsep_test_v2(
                  i,
                  test_eq,
                  current_monitor,
                  engine = engine,
                  nimble_samplers = nimble_samplers,
                  quiet = quiet,
                  original_data = original_data,
                  hierarchical_info = hierarchical_info,
                  random_terms = random_terms,
                  equations = equations,
                  family = family,
                  structure = structure,
                  levels = levels,
                  hierarchy = hierarchy,
                  link_vars = link_vars,
                  fix_residual_variance = fix_residual_variance,
                  latent_method = latent_method,
                  n.chains = n.chains,
                  n.iter = n.iter,
                  n.burnin = n.burnin,
                  n.thin = n.thin,
                  n.adapt = n.adapt,
                  ic_recompile = ic_recompile,
                  latent = latent,
                  random = random,
                  id_col = id_col,
                  variability = variability
                )
              }, error = function(e) {
                # Return error info to parent session for reporting
                list(
                  error = conditionMessage(e),
                  test_idx = i,
                  equation = paste(deparse(test_eq), collapse = " ")
                )
              })
            }
          )
        }

        # Merge parallel results into new_results_list
        for (j in seq_along(tests_to_run_indices)) {
          idx <- tests_to_run_indices[j]
          res <- par_results[[j]]
          
          # Check if the result is an error object from tryCatch
          if (is.list(res) && !is.null(res$error)) {
            if (!quiet) {
              warning(sprintf(
                "D-sep test %d/%d skipped due to error: %s\n  Equation: %s",
                res$test_idx, length(dsep_tests), res$error, res$equation
              ))
            }
            new_results_list[[idx]] <- NULL
          } else {
            new_results_list[[idx]] <- res
          }
        }
      } else {
        # Sequential execution
        # Loop over only the needed indices
        seq_results <- lapply(tests_to_run_indices, function(i) {
          test_eq <- dsep_tests[[i]]
          # Monitor betas for ALL predictors in the d-sep equation
          current_monitor <- "interpretable"

          if (!quiet) {
            message(sprintf(
              "  Test %d/%d: %s",
              i,
              length(dsep_tests),
              deparse(test_eq)
            ))
          }
          new_results_list[[i]] <- tryCatch(
            run_single_dsep_test_v2(
              i,
              test_eq,
              current_monitor,
              engine = engine,
              nimble_samplers = nimble_samplers,
              quiet = quiet,
              original_data = original_data,
              hierarchical_info = hierarchical_info,
              random_terms = random_terms,
              equations = equations,
              family = family,
              structure = structure,
              levels = levels,
              hierarchy = hierarchy,
              link_vars = link_vars,
              fix_residual_variance = fix_residual_variance,
              latent = latent,
              latent_method = latent_method,
              n.chains = n.chains,
              n.iter = n.iter,
              n.burnin = n.burnin,
              n.thin = n.thin,
              n.adapt = n.adapt,
              ic_recompile = ic_recompile,
              random = random,
              id_col = id_col,
              variability = variability
            ),
            error = function(e) {
              if (!quiet) {
                warning(sprintf(
                  "D-sep test %d/%d skipped due to error: %s\n  Equation: %s",
                  i, length(dsep_tests), conditionMessage(e), deparse(test_eq)
                ))
              }
              NULL  # Return NULL for this test
            }
          )
        })

        # Assign sequential results to the main list
        for (j in seq_along(tests_to_run_indices)) {
          idx <- tests_to_run_indices[j]
          new_results_list[[idx]] <- seq_results[[j]]
        }
      }
    } # End if tests_to_run > 0

    results <- new_results_list

    # Combine results
    combined_models <- list()

    for (res_item in results) {
      if (is.null(res_item)) {
        next
      }

      tryCatch(
        {
          samples <- res_item$samples
          param_map <- res_item$param_map
          model_string <- res_item$model
          i <- res_item$test_index

          # Store model for this test
          combined_models[[i]] <- model_string

          # Rename parameters to include equation index to avoid collisions
          # e.g., betaRS becomes betaRS_1 for equation 1, betaRS_2 for equation 2
          for (ch in seq_along(samples)) {
            chain <- samples[[ch]]
            colnames_orig <- colnames(chain)

            # Add suffix _i to all beta, alpha, lambda, tau, rho parameters
            new_colnames <- sapply(colnames_orig, function(name) {
              if (grepl("^(beta|alpha|lambda|tau|rho|sigma)", name)) {
                paste0(name, "_", i)
              } else {
                name
              }
            })

            colnames(chain) <- new_colnames
            samples[[ch]] <- chain
          }

          # Update parameter_map to reflect new names
          if (!is.null(param_map) && nrow(param_map) > 0) {
            param_map$parameter <- paste0(param_map$parameter, "_", i)
          }

          # Combine samples (cbind chains)
          if (is.null(combined_samples)) {
            combined_samples <- samples
          } else {
            # Check if dimensions match
            if (coda::niter(combined_samples) != coda::niter(samples)) {
              stop("MCMC iteration mismatch between d-sep tests")
            }
            # Combine chains: for each chain, cbind the variables
            new_samples <- coda::mcmc.list()
            for (ch in 1:coda::nchain(combined_samples)) {
              # Combine matrices
              mat1 <- combined_samples[[ch]]
              mat2 <- samples[[ch]]
              # All columns from mat2 should be new (due to renaming)
              new_mat <- cbind(mat1, mat2)
              new_samples[[ch]] <- coda::mcmc(
                new_mat,
                start = stats::start(mat1),
                thin = coda::thin(mat1)
              )
            }
            combined_samples <- new_samples
          }

          # Combine parameter maps
          if (is.null(combined_map)) {
            combined_map <- param_map
          } else {
            combined_map <- rbind(combined_map, param_map)
          }
        },
        error = function(e) {
          if (!quiet) {
            warning(paste(
              "D-sep result combination error for test index",
              res_item$test_index,
              ":",
              e$message
            ))
          }
        }
      )
    }

    # Return combined result
    result <- list(
      samples = combined_samples,
      parameter_map = combined_map,
      models = combined_models,
      dsep = TRUE,
      dsep_tests = dsep_tests,
      dsep_results = results, # Store individual test results for summary
      induced_correlations = induced_cors,
      equations = equations,
      data = data,
      original_data = original_data,
      family = family,
      categorical_vars = attr(data, "categorical_vars"),
      poly_terms = all_poly_terms
    )
    class(result) <- "because"
    return(result)
  }

  # Handle latent variable method
  if (!is.null(latent)) {
    latent_method <- match.arg(latent_method, c("correlations", "explicit"))

    # Force MAG approach when doing d-separation testing
    if (dsep && latent_method == "explicit") {
      if (!quiet) {
        message(
          "Note: d-separation testing with latent variables requires MAG approach. ",
          "Using latent_method = 'correlations'."
        )
      }
      latent_method <- "correlations"
    }

    if (latent_method == "correlations") {
      # MAG approach: marginalize latents, use induced correlations
      # If not already computed by dsep, compute now
      if (is.null(induced_cors)) {
        dsep_result <- because_dsep(
          equations,
          latent = latent,
          random_terms = random_terms,
          hierarchical_info = hierarchical_info,
          family = family,
          quiet = !dsep
        )
        induced_cors <- dsep_result$correlations
      }

      # Filter out equations involving latent variables
      # Handle equations involving latent variables
      original_eq_count <- length(equations)
      equations <- lapply(equations, function(eq) {
        vars <- all.vars(eq)
        if (any(vars %in% latent)) {
          # Remove latent variables from the formula
          # We construct a new formula string excluding latent vars
          rhs <- labels(terms(eq))
          keep_terms <- rhs[!rhs %in% latent]

          if (length(keep_terms) == 0) {
            # Becomes intercept-only model
            new_eq <- as.formula(paste(as.character(eq)[2], "~ 1"))
          } else {
            # Keep observed predictors
            new_eq <- as.formula(paste(
              as.character(eq)[2],
              "~",
              paste(keep_terms, collapse = " + ")
            ))
          }
          return(new_eq)
        }
        return(eq)
      })

      # We no longer filter them out, so we don't report "removed equations"
      if (!quiet) {
        message(
          "Using MAG approach: marginalized latent variables from structural equations."
        )
      }

      # Check if variables with induced correlations need intercept-only models
      # This happens when a variable is involved in an induced correlation
      # but is not a response variable in any remaining equation
      if (length(induced_cors) > 0) {
        # Get all variables from induced correlations
        vars_with_correlations <- unique(unlist(induced_cors))

        # Get response variables from remaining equations
        response_vars <- sapply(equations, function(eq) {
          as.character(eq)[2]
        })

        # Find variables that need intercept models
        vars_needing_intercept <- setdiff(
          vars_with_correlations,
          response_vars
        )

        if (length(vars_needing_intercept) > 0) {
          # Create intercept-only models: X ~ 1
          intercept_equations <- lapply(vars_needing_intercept, function(v) {
            as.formula(paste(v, "~ 1"))
          })

          # Add to equations list
          equations <- c(equations, intercept_equations)

          if (!quiet) {
            message(
              "Created intercept-only models for ",
              length(vars_needing_intercept),
              " variable(s) with induced correlations: ",
              paste(vars_needing_intercept, collapse = ", ")
            )
          }
        }
      }

      if (!quiet && length(induced_cors) > 0) {
        message(
          "Estimating ",
          length(induced_cors),
          " induced correlation(s) from latent variable(s)"
        )
      }
    } else {
      # Explicit approach: keep all equations, don't use induced correlations
      induced_cors <- NULL

      if (!quiet) {
        message("Using explicit latent variable modeling")
      }
    }
  }

  # Auto-expand categorical variables in equations
  if (!is.null(attr(data, "categorical_vars"))) {
    categorical_vars <- attr(data, "categorical_vars")
    new_equations <- list()

    # Use for loop instead of lapply to safely modify data and collect equations
    for (idx in seq_along(equations)) {
      eq <- equations[[idx]]
      # Parse formula to get all variables
      vars <- all.vars(eq)

      # Check if any predictors are categorical
      for (var in vars) {
        if (var %in% names(categorical_vars)) {
          # Check if var is the response (LHS)
          lhs_var <- all.vars(eq[[2]])
          if (var %in% lhs_var) {
            # Skip expansion if it's the response
            next
          }

          # Get dummy variable names
          levels <- categorical_vars[[var]]$levels
          dummies <- categorical_vars[[var]]$dummies

          # If the parent variable is being imputed (is in response_vars_with_na),
          # we MUST define the dummies deterministically in JAGS to link them.
          if (var %in% response_vars_with_na) {
            for (k in 2:length(levels)) {
              dummy_name <- paste0(var, "_", levels[k])

              # Create deterministic equation: dummy ~ I(var == k)
              det_eq_str <- sprintf("%s ~ I(%s == %d)", dummy_name, var, k)
              det_eq <- stats::as.formula(det_eq_str)

              # Only add if not already present
              eq_exists <- any(sapply(c(equations, new_equations), function(e) {
                deparse(e) == deparse(det_eq)
              }))

              if (!eq_exists) {
                new_equations <- c(new_equations, list(det_eq))
                # Also remove the dummy from the data list so JAGS uses the definition
                if (dummy_name %in% names(data)) {
                  data[[dummy_name]] <- NULL
                }
              }
            }
          }

          # [FIX] Skip substitution if we are inside a deterministic identity definition
          # for an imputed factor (e.g., Rate_2 ~ I(Rate == 2)).
          # This prevents circular expansion of 'Rate' into its own dummies.
          is_identity_mapping <- length(vars) == 2 && grepl(sprintf("I\\(%s == [0-9]+\\)", var), deparse(eq))
          if (is_identity_mapping) {
            next
          }

          # Convert formula to character for manipulation
          eq_str <- paste(deparse(eq), collapse = " ")

          # Replace categorical variable with its dummies (wrapped in parentheses for proper interaction expansion like A*B)
          pattern <- paste0("\\b", var, "\\b")
          replacement <- paste0("(", paste(dummies, collapse = " + "), ")")
          eq_str <- gsub(pattern, replacement, eq_str)

          # Convert back to formula
          eq <- as.formula(eq_str)
          equations[[idx]] <- eq

          if (!quiet) {
            message(sprintf(
              "Expanded '%s' to: %s",
              var,
              paste(dummies, collapse = ", ")
            ))
          }
        }
      }
    }
    # Add the new deterministic equations
    equations <- c(equations, new_equations)
  }

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
    family = family,
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
    # We call because_dsep just for its side effect (printing MAG info).
    # We wrap it in tryCatch to ensure it doesn't block the main run if it fails.
    tryCatch(
      {
        if (length(equations) > 0) {
          invisible(because_dsep(
            equations,
            latent = latent,
            random_terms = random_terms,
            quiet = quiet
          ))
        }
      },
      error = function(e) {
        # Silently skip if MAG display fails - not critical for main run
        if (!quiet) {
          message("(MAG structure display skipped)")
        }
      }
    )
    message("---------------------------------------")
  }

  # Monitor parameters
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
          "(?:logit|log|cloglog|probit)?\\(?\\s*(\\w+)(?:\\[.*\\])?\\)?\\s*(?:<-|~)",
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
      "\\b([a-zA-Z0-9_]+)\\s*\\[1:N\\]\\s*~",
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

        m_obj <- nimble::nimbleModel(
          code = nimble_code,
          constants = data,
          inits = nimble_inits
        )
        
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
        m_obj
      },
      error = function(e) {
        if (!quiet) {
          message("\nCRITICAL NIMBLE ERROR during model initialization:")
          message(e$message)
        }
        stop(e)
      }
    )

    # Configure MCMC
    mcmc_conf <- nimble::configureMCMC(
      nimble_model,
      monitors = monitor,
      enableWAIC = WAIC
    )

    # Harden NIMBLE sampler assignments.
    # Called via the top-level exported function so the code is always
    # taken from the INSTALLED package, not from this closure.
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

        worker_model <- nimble::nimbleModel(
          code = nimble_code,
          constants = data,
          inits = curr_inits
        )

        worker_conf <- nimble::configureMCMC(
          worker_model,
          monitors = monitor,
          enableWAIC = WAIC
        )
        
        # Harden worker MCMC samplers via the installed package function.
        # CRITICAL: calling by package::function rather than using the closure
        # ensures workers always run the CURRENT installed version of because,
        # not the stale closure that was exported from the parent session.
        because::nimble_harden_samplers(
          worker_conf,
          family          = family,
          nimble_samplers = nimble_samplers,
          quiet           = TRUE
        )

        worker_mcmc <- nimble::buildMCMC(worker_conf)
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
          family = family,
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

      # Format into mcmc.list
      samples <- coda::mcmc.list(lapply(chain_results, function(x) x))
    } else {
      # Sequential execution (Existing logic)
      if (!quiet) {
        message(sprintf(
          "Sampling %d chains sequentially via NIMBLE...",
          n.chains
        ))
      }

      nimble_samples <- nimble::runMCMC(
        compiled_mcmc,
        niter = n.iter,
        nburnin = n.burnin,
        nchains = n.chains,
        thin = n.thin,
        samplesAsCodaMCMC = TRUE,
        summary = FALSE
      )

      if (n.chains == 1) {
        samples <- coda::mcmc.list(nimble_samples)
      } else {
        samples <- coda::mcmc.list(lapply(nimble_samples, function(x) x))
      }
    }
    model <- nimble_model # Store the Rmodel object
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

        model <- rjags::jags.model(
          model_file,
          data = data,
          inits = inits_list,
          n.chains = 1,
          n.adapt = n.adapt,
          quiet = quiet
        )

        # Burn-in
        if (n.burnin > 0) {
          update(model, n.iter = n.burnin)
        }

        # Sample
        if (n.iter > n.burnin) {
          samples <- rjags::coda.samples(
            model,
            variable.names = monitor,
            n.iter = n.iter - n.burnin,
            thin = n.thin
          )
        } else {
          samples <- NULL
        }

        return(list(samples = samples, model = model))
      }

      # Export necessary objects to cluster
      parallel::clusterExport(cl, c("run_single_chain"), envir = environment())

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

      model <- tryCatch(
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
          stop(e)
        }
      )
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

      # Sample posterior
      if (n.iter > n.burnin) {
        samples <- rjags::coda.samples(
          model,
          variable.names = monitor,
          n.iter = n.iter - n.burnin,
          thin = n.thin
        )
      } else {
        samples <- NULL
      }
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
      family = family,
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
    stacked_data = if (exists("stack_res") && stack_res$is_stacked) {
      stack_res$data
    } else {
      NULL
    }
  )

  # Combine with basic model info
  result$model <- model_string
  result$samples <- samples
  result$parameter_map <- parameter_map
  result$data <- data
  result$original_data <- original_data
  result$family <- family
  result$categorical_vars <- attr(data, "categorical_vars")
  result$poly_terms <- all_poly_terms
  result$equations <- equations

  # --- Result Enrichment ---
  # If we have a structure, try to extract labels for ordering
  if (!is.null(structure)) {
    result$species_order <- get_order_labels_hook(structure)
  } else if (!is.null(id_col) && is.data.frame(original_data)) {
    # If no structure but ID col provided
    result$species_order <- as.character(original_data[[id_col]])
  }

  # Assign class immediately (needed for print/summary/waic methods)
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
    result$WAIC <- because_waic(result)
  }

  # Preserve hierarchical metadata for diagnostics even if data was flat
  result$hierarchical_info <- hierarchical_info
  
  return(result)
}

#' @noRd
run_single_dsep_test_v2 <- function(
  i,
  test_eq,
  monitor_params,
  engine = "jags",
  nimble_samplers = NULL,
  quiet = FALSE,
  original_data = NULL,
  hierarchical_info = NULL,
  random_terms = list(),
  equations = list(),
  family = NULL,
  structure = NULL,
  levels = NULL,
  hierarchy = NULL,
  link_vars = NULL,
  fix_residual_variance = NULL,
  latent = NULL,
  latent_method = "correlations",
  n.chains = 3,
  n.iter = 12500,
  n.burnin = 2500,
  n.thin = 10,
  n.adapt = 2500,
  ic_recompile = FALSE,
  random = NULL,
  id_col = NULL,
  variability = NULL
) {
  if (!quiet) {
    message(paste("D-sep test eq:", deparse(test_eq)))
  }

  # Select appropriate dataset for this test
  test_data <- original_data

  if (!is.null(hierarchical_info)) {
    # Extract variables from this test equation
    test_vars <- all.vars(test_eq)
    
    # [FIX 2026-04-13] we NO LONGER strip link_vars from test_vars.
    # infer_variable_level() now has a fallback to search data columns,
    # and get_data_for_variables() needs these IDs to correctly assemble/join levels.
    
    # [FIX] Add random effect grouping variables to test_vars
    # Otherwise get_data_for_variables removes them, causing "Unknown variable N_SiteID"
    # [FIX] categorical dummy variables
    # We need their PARENT variables to fetch the data, then recreate dummies manually.
    dummies_to_create <- list()
    cat_vars <- NULL

    if (!is.null(attr(original_data, "categorical_vars"))) {
      cat_vars <- attr(original_data, "categorical_vars")

      # Check all cat vars to see if their dummies are needed (or if parent is needed)
      current_vars <- test_vars

      for (parent_var in names(cat_vars)) {
        dummies <- cat_vars[[parent_var]]$dummies

        # If parent variable is in the test variables, we need to ensure we can recreate expected dummies
        if (parent_var %in% current_vars) {
          dummies_to_create[[parent_var]] <- dummies
        }
      }
    }

    # Get appropriate dataset for these variables (dummies NOT included in request)
    test_data <- get_data_for_variables(
      test_vars,
      original_data,
      hierarchical_info$levels,
      hierarchical_info$hierarchy,
      hierarchical_info$link_vars,
      equations = equations,
      latent = latent
    )
    
    # [FIX] Preserve categorical metadata for K detection in sub-models
    if (!is.null(cat_vars)) {
      attr(test_data, "categorical_vars") <- cat_vars
    }

    # [FIX] Recreate dummy variables in test_data
    for (parent_var in names(dummies_to_create)) {
      if (parent_var %in% names(test_data)) {
        dummies <- dummies_to_create[[parent_var]]
        vals <- test_data[[parent_var]]

        # Recreate dummies: sex_m = as.integer(sex == 2) etc.
        cat_vars_attr <- attr(original_data, "categorical_vars")
        levels_map <- if (
          !is.null(cat_vars_attr) && parent_var %in% names(cat_vars_attr)
        ) {
          cat_vars_attr[[parent_var]]$levels
        } else {
          NULL
        }
        if (is.null(levels_map)) {
          next
        }

        # Create all dummies for this parent
        type <- cat_vars[[parent_var]]$type
        
        if (!is.null(type) && type == "ordered") {
          c_mat <- cat_vars[[parent_var]]$contrasts
          vals <- test_data[[parent_var]]
          
          # Handle values which might be numeric indices or strings
          match_idx <- if (is.numeric(vals)) {
            vals
          } else {
            match(vals, levels_map)
          }
          
          if (any(is.na(match_idx))) {
            print("NA found in match_idx!")
            print(head(vals))
            print(levels_map)
          }

          for (k in seq_along(dummies)) {
            expected_dummy <- dummies[k]
            if (is.null(test_data[[expected_dummy]]) || all(is.na(test_data[[expected_dummy]]))) {
              test_data[[expected_dummy]] <- c_mat[match_idx, k]
            }
          }
        } else {
          for (k in 2:length(levels_map)) {
            expected_dummy <- dummies[k - 1]
  
            # Check if it was extracted properly. If not, recreate it!
            # We only recreate if it's missing or NA
            if (
              is.null(test_data[[expected_dummy]]) ||
                all(is.na(test_data[[expected_dummy]]))
            ) {
              vals <- test_data[[parent_var]]
  
              # Robust check matching preprocess_categorical_vars logic
              is_match <- if (is.numeric(vals)) {
                vals == k
              } else {
                vals == levels_map[k]
              }
  
              test_data[[expected_dummy]] <- as.integer(is_match)
            }
          }
        }
      }
    }

    # [FIX] Restore categorical_vars attribute dropped by merge/get_data
    if (!is.null(attr(original_data, "categorical_vars"))) {
      attr(test_data, "categorical_vars") <- attr(
        original_data,
        "categorical_vars"
      )
    }

    if (!quiet) {
      message(
        "  Hierarchical data: using ",
        nrow(test_data),
        " observations for this test"
      )
    }
  }

  # [REVISED] Populate dsep_equations FIRST, then filter sub_family/sub_variability
  dsep_equations <- list(test_eq)

  test_eq_vars <- all.vars(test_eq)
  # [NEW 2025-12-21] For occupancy models, we must include supporting equations
  # (detection models and models for latent predictors)
  # [Fixed duplicate line]

  # Extension Hook: Expand d-separation equations (e.g. detection models)
  dsep_equations <- dsep_equations_hook(
    family,
    equations,
    dsep_equations,
    test_eq = test_eq
  )

  # Now filter variability/family based on ALL variables in dsep_equations
  # [Manual Fix: variability needs to be handled if it's there, but here we only have family]
  all_dsep_vars <- unique(unlist(lapply(dsep_equations, all.vars)))
  all_dsep_responses <- unique(unlist(lapply(dsep_equations, function(eq) as.character(eq)[2])))
  
  all_dsep_vars_clean <- unique(c(
    all_dsep_vars,
    sub("^p_", "", all_dsep_vars),
    sub("^psi_", "", all_dsep_vars),
    sub("^z_", "", all_dsep_vars)
  ))
  
  all_dsep_responses_clean <- unique(c(
    all_dsep_responses,
    sub("^p_", "", all_dsep_responses),
    sub("^psi_", "", all_dsep_responses),
    sub("^z_", "", all_dsep_responses)
  ))

  # Note: sub_variability logic removed for brevity if not strictly needed in this context
  # but we'll try to extract what we can from arguments
  sub_family <- if (!is.null(family)) {
    family[names(family) %in% all_dsep_responses_clean]
  } else {
    NULL
  }

  # [User Request] Auto-detect binomial for binary response
  test_resp <- as.character(test_eq)[2]
  if (!is.null(test_data[[test_resp]])) {
    # Try to find target column in test_data (which might be a list or df)
    vals <- if (is.list(test_data) && !is.data.frame(test_data)) {
      # find which element contains it
      found_vals <- NULL
      for (lvl in names(test_data)) {
        if (test_resp %in% names(test_data[[lvl]])) {
          found_vals <- test_data[[lvl]][[test_resp]]
          break
        }
      }
      found_vals
    } else {
      test_data[[test_resp]]
    }

    if (!is.null(vals)) {
      u_vals <- unique(na.omit(vals))
      if (length(u_vals) <= 2 && all(u_vals %in% c(0, 1))) {
        if (is.null(sub_family)) {
          sub_family <- list()
        }
        if (is.na(sub_family[test_resp])) {
          sub_family[[test_resp]] <- "binomial"
        }
      }
    }
  }

  if (length(sub_family) == 0) {
    sub_family <- NULL
  }

  # Choose what data to pass to the dsep sub-fit:
  # [FIX] Never pass the flattened joint table (hierarchical_info$data) to sub-models!
  # This avoids inheriting the observation-level resolution and prevents inflation.
  # We pass original_data (the raw list) to trigger resolution-locked assembly.
  dsep_data_to_pass <- if (!is.null(original_data)) {
    original_data
  } else {
    test_data
  }

  # Extension Hook: Filter structure
  sub_structure <- dsep_tree_hook(structure, test_eq, hierarchical_info, levels)

  # Call because recursively
  # Use do.call and filtering to handle potential version conflicts on worker nodes
  bec_args <- names(formals(because))
  call_args <- list(
    data = dsep_data_to_pass,
    structure = sub_structure,
    equations = dsep_equations,
    monitor = monitor_params,
    n.chains = n.chains,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = n.thin,
    DIC = FALSE,
    WAIC = FALSE,
    n.adapt = n.adapt,
    quiet = TRUE,
    dsep = FALSE,
    family = sub_family,
    fix_residual_variance = fix_residual_variance,
    latent = latent,
    latent_method = latent_method,
    parallel = FALSE,
    n.cores = 1,
    cl = NULL,
    ic_recompile = ic_recompile,
    random = random,
    levels = levels,
    hierarchy = hierarchy,
    link_vars = link_vars,
    structure_multi = hierarchical_info$structure_multi,
    structure_levels = hierarchical_info$structure_levels,
    id_col = id_col,
    variability = variability
  )

  # Only add NIMBLE arguments if the version of because() on this node supports them
  if ("engine" %in% bec_args) {
    call_args$engine <- engine
    call_args$nimble_samplers <- nimble_samplers
  }

  fit <- do.call(because, call_args)

  # Extract samples, map, and model
  samples <- fit$samples
  model_string <- fit$model
  param_map <- fit$parameter_map

  # Update equation index in parameter map to match the d-sep test index
  param_map$equation_index <- i

  list(
    samples = samples,
    param_map = param_map,
    model = model_string,
    test_index = i
  )
}

#' Preprocess categorical variables (character/factor) to integer codes and dummies
#'
#' @param data A data.frame or list of data.frames
#' @param quiet Logical; whether to suppress informational messages
#' @return The modified data object with categorical_vars attribute
#' @keywords internal
preprocess_categorical_vars <- function(
  data,
  target_vars = NULL,
  dummy_vars = NULL,
  exclude_cols = NULL,
  quiet = FALSE,
  expand_ordered = FALSE
) {
  if (is.null(data)) {
    return(NULL)
  }

  # Recursively process lists of data frames (hierarchical data)
  if (is.list(data) && !is.data.frame(data)) {
    all_cat_vars <- list()

    # Process each element (usually levels in hierarchy)
    for (i in seq_along(data)) {
      if (is.data.frame(data[[i]])) {
        processed <- preprocess_categorical_vars(
          data[[i]],
          target_vars = target_vars,
          dummy_vars = dummy_vars,
          exclude_cols = exclude_cols,
          quiet = TRUE, # Suppress output for recursive calls to avoid noise
          expand_ordered = expand_ordered
        )
        data[[i]] <- processed
        # Collect categorical vars metadata
        level_cat_vars <- attr(processed, "categorical_vars")
        if (!is.null(level_cat_vars)) {
          # Use utils::modifyList if available, or manual merge
          for (name in names(level_cat_vars)) {
            all_cat_vars[[name]] <- level_cat_vars[[name]]
          }
        }
      }
    }

    # Attach merged metadata to the top-level list
    attr(data, "categorical_vars") <- all_cat_vars
    return(data)
  }

  # Process single data frame
  if (!is.data.frame(data)) {
    return(data)
  }

  # Determine which columns to check
  check_cols <- if (!is.null(target_vars)) {
    intersect(names(data), target_vars)
  } else {
    names(data)
  }

  if (length(check_cols) == 0) {
    return(data)
  }

  # Detect categorical variables: either by type or if they already have metadata
  cat_metadata <- attr(data, "categorical_vars")
  char_cols <- sapply(check_cols, function(col) {
    x <- data[[col]]
    is.character(x) || is.factor(x) || (!is.null(cat_metadata) && col %in% names(cat_metadata))
  })

  if (any(char_cols)) {
    categorical_vars <- list()
    if (!is.null(attr(data, "categorical_vars"))) {
      categorical_vars <- attr(data, "categorical_vars")
    }

    col_names <- check_cols[char_cols]
    for (col in col_names) {
      if (!is.null(exclude_cols) && col %in% exclude_cols) {
        next
      }

      # [FIX] Preserve existing categorical metadata if present (to keep levels consistent in subsets)
      existing_metadata <- categorical_vars[[col]]
      if (!is.null(existing_metadata)) {
        levels <- existing_metadata$levels
        is_ord_prev <- !is.null(existing_metadata$type) && existing_metadata$type == "ordered"
        
        # If already integer/numeric, it's likely already encoded from a previous call.
        # Use integer-to-label mapping to avoid wiping the data with NAs.
        if (is.numeric(data[[col]])) {
          f_vals <- factor(data[[col]], levels = seq_along(levels), labels = levels, ordered = is_ord_prev)
        } else {
          # Use existing levels to ensure integer encoding (1, 2, 3...) remains consistent
          # during d-separation tests on data subsets.
          f_vals <- factor(data[[col]], levels = levels, ordered = is_ord_prev)
        }
      } else {
        # New variable: Convert to factor to discover levels
        f_vals <- factor(data[[col]])
        levels <- levels(f_vals)
      }

      if (length(levels) < 2) {
        if (!quiet) {
          warning(sprintf(
            "Variable '%s' has < 2 levels. Converting to numeric constant.",
            col
          ))
        }
        data[[col]] <- as.numeric(f_vals)
      } else {
        # Store metadata for model expansion
        is_ord <- is.ordered(f_vals)

        if (is.null(existing_metadata)) {
          if (is_ord) {
            if (expand_ordered) {
              c_mat <- stats::contr.poly(length(levels))
              if (ncol(c_mat) > 2) {
                c_mat <- c_mat[, 1:2, drop = FALSE]
              }
              c_names <- colnames(c_mat)
              c_names[c_names == ".L"] <- "L"
              c_names[c_names == ".Q"] <- "Q"
              c_names[c_names == ".C"] <- "C"
              c_names <- gsub("\\^", "pow", c_names)

              categorical_vars[[col]] <- list(
                levels = levels,
                reference = "Polynomial Contrast",
                dummies = paste0(col, "_", c_names),
                type = "ordered",
                contrasts = c_mat
              )
            } else {
              # Default: Linear only (centered integers)
              categorical_vars[[col]] <- list(
                levels = levels,
                reference = "Numeric (Centered)",
                dummies = paste0(col, "_L"),
                type = "ordered"
              )
            }
          } else {
            categorical_vars[[col]] <- list(
              levels = levels,
              reference = levels[1],
              dummies = paste0(col, "_", levels[-1]),
              type = "unordered"
            )
          }
        } else {
          # Update existing metadata if needed
          if (is_ord && is.null(existing_metadata$contrasts) && expand_ordered) {
            c_mat <- stats::contr.poly(length(levels))
            if (ncol(c_mat) > 2) {
              c_mat <- c_mat[, 1:2, drop = FALSE]
            }
            c_names <- colnames(c_mat)
            c_names[c_names == ".L"] <- "L"
            c_names[c_names == ".Q"] <- "Q"
            c_names[c_names == ".C"] <- "C"
            c_names <- gsub("\\^", "pow", c_names)

            categorical_vars[[col]]$reference <- "Polynomial Contrast"
            categorical_vars[[col]]$dummies <- paste0(col, "_", c_names)
            categorical_vars[[col]]$type <- "ordered"
            categorical_vars[[col]]$contrasts <- c_mat
          } else if (is_ord && is.null(existing_metadata$dummies) && !expand_ordered) {
            categorical_vars[[col]]$reference <- "Numeric (Centered)"
            categorical_vars[[col]]$dummies <- paste0(col, "_L")
            categorical_vars[[col]]$type <- "ordered"
          }
        }

        # Convert to integer codes for JAGS
        data[[col]] <- as.integer(f_vals)

        # Generate Dummy Variables explicitly ONLY if requested
        if (is.null(dummy_vars) || col %in% dummy_vars) {
          dummies <- categorical_vars[[col]]$dummies

          if (!quiet && length(levels) > 500) {
            message(sprintf(
              "Generating %d dummy variables for '%s'... this may take a moment.",
              length(levels) - 1,
              col
            ))
          }
          if (is_ord) {
            if (expand_ordered) {
              c_mat <- categorical_vars[[col]]$contrasts
              for (k in seq_along(dummies)) {
                data[[dummies[k]]] <- c_mat[data[[col]], k]
              }
            } else {
              # Linear centered (dummy name is col_L)
              raw_codes <- data[[col]]
              data[[paste0(col, "_L")]] <- raw_codes - mean(raw_codes, na.rm = TRUE)
            }
          } else {
            for (k in 2:length(levels)) {
              lev_name <- levels[k]
              dummy_col_name <- dummies[k - 1]
              # Create binary column: 1 if matches this level, 0 if not (NAs remain NA)
              data[[dummy_col_name]] <- as.numeric(data[[col]] == k)
            }
          }
        }

        if (!quiet) {
          message(sprintf(
            "Converted categorical '%s' to integers (1..%d). Reference: '%s'%s",
            col,
            length(levels),
            levels[1],
            if (is.null(dummy_vars) || col %in% dummy_vars) {
              ""
            } else {
              " (Dummies skipped)"
            }
          ))
        }
      }
    }
    attr(data, "categorical_vars") <- categorical_vars
  }

  return(data)
}
