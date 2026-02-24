#' @importFrom stats dist
#' @importFrom utils read.csv
#' @importFrom ape vcv
#' @importFrom methods is
#' @importFrom rjags jags.model
#'
#' @title Run a Bayesian Structural Equation Model
#'
#' @description
#' This function fits a Bayesian ...
#' @param equations A list of model formulas describing the structural equation model.
#' @param id_col Character string specifying the column name in a data.frame containing
#'   unit identifiers (species, individuals, sites, etc.). This is used to:
#'   \itemize{
#'     \item Match data rows to tree tip labels (for phylogenetic models)
#'     \item Link data to external spatial or custom covariance matrices.
#'   }
#'   **Note**: For standard random effects models (e.g. \code{random = ~(1|species)}) where no external structure
#'   (like a tree) is provided, this argument is **not required**. The grouping column is read directly from the data.
#'
#'   If \code{NULL} (default): uses meaningful row names if available.
#'   Ignored when \code{data} is already a list.
#' @param structure The covariance structure for the model. Accepts:
#'   \itemize{
#'     \item \code{"phylo"} object: Phylogenetic tree (Standard PGLS/PhyloSEM).
#'     \item \code{"multiPhylo"} object: List of trees (incorporates phylogenetic uncertainty).
#'     \item \code{NULL}: Independent model (Standard SEM, no covariance structure).
#'     \item \code{matrix}: Custom covariance or precision matrix (e.g., spatial connectivity, kinship).
#'   }
#' @param tree (Deprecated alias for \code{structure}). A single phylogenetic tree of class
#'   \code{"phylo"} or a list of trees. Use \code{structure} instead for new code.
#' @param monitor Parameter monitoring mode. Options:
#'   \itemize{
#'     \item \code{"interpretable"} (default): Monitor only scientifically meaningful parameters:
#'           intercepts (alpha), regression coefficients (beta), phylogenetic signals (lambda) for
#'          responses, and WAIC terms. Excludes variance components (tau) and auxiliary predictor parameters.
#'     \item \code{"all"}: Monitor all model parameters including variance components and implicit equation parameters.
#'     \item Custom vector: Provide a character vector of specific parameter names to monitor.
#'     \item \code{NULL}: Auto-detect based on model structure (equivalent to "interpretable").
#'   }
#' @param n.chains Number of MCMC chains (default = 3).
#' @param n.iter Total number of MCMC iterations (default = 12500).
#' @param n.burnin Number of burn-in iterations (default = n.iter / 5).
#' @param n.thin Thinning rate (default = 10).
#' @param DIC Logical; whether to compute DIC using \code{dic.samples()} (default = TRUE).
#'   **Note**: DIC penalty will be inflated for models with measurement error or repeated measures
#'   because latent variables are counted as parameters (penalty ~ structural parameters + N).
#'   For model comparison, use WAIC or compare mean deviance across models with similar structure.
#' @param WAIC Logical; whether to sample values for WAIC and deviance (default = FALSE).
#'   WAIC is generally more appropriate than DIC for hierarchical models with latent variables.
#' @param n.adapt Number of adaptation iterations (default = n.iter / 5).
#' @param quiet Logical; suppress JAGS output (default = FALSE).
#' @param verbose Logical; if \code{TRUE}, print generated JAGS model code and data names (default = FALSE).
#' @param dsep Logical; if \code{TRUE}, monitor only the first beta in each structural equation (used for d-separation testing).
#' @param variability Optional specification for variables with measurement error or within-species variability.
#'   \strong{Global Setting}:
#'   \itemize{
#'     \item \code{"reps"}: Applies repeat-measures modeling to \strong{all} continuous variables in the equations (except grouping variables). Expects \code{X_obs} matrix or long-format data.
#'     \item \code{"se"}: Applies measurement error modeling to \strong{all} continuous variables. Expects \code{X_se} columns.
#'   }
#'
#'   \strong{Manual Specification} (Named Vector/List):
#'   \itemize{
#'     \item Simple: \code{c(X = "se", Y = "reps")} - mixed types
#'     \item Custom columns: \code{list(X = list(type = "se", se_col = "X_sd"))}
#'     \item For SE: \code{se_col} (SE column), \code{mean_col} (mean column, optional)
#'     \item For reps: \code{obs_col} (observations matrix column)
#'   }
#'
#'   \strong{Auto-Detection}:
#'   If not specified, the package attempts to detect variability based on column names:
#'   \itemize{
#'     \item \code{X_se} -> type="se"
#'     \item \code{X_obs} or matrix column -> type="reps"
#'   }
#' @param family Optional named character vector specifying the family/distribution for response variables.
#'   Default is "gaussian" for all variables. Supported values:
#'   \itemize{
#'     \item "gaussian" (default)
#'     \item "binomial" (binary data)
#'     \item "multinomial" (unordered categorical > 2 levels)
#'     \item "ordinal" (ordered categorical > 2 levels)
#'     \item "poisson" (count data)
#'     \item "negbinomial" (overdispersed count data)
#'     \item "zip" (zero-inflated poisson): Models excess zeros with probability \code{psi} and counts with mean \code{lambda}.
#'     \item "zinb" (zero-inflated negative binomial): Models excess zeros with probability \code{psi} and overdispersed counts with mean \code{mu} and size \code{r}.
#'     \item "occupancy": Single-season site-occupancy model.
#'           State process: \code{z ~ Bernoulli(psi)}.
#'           Observation process: \code{y ~ Bernoulli(z * p)}.
#'           Requires data to be a detection history matrix (sites x visits).
#'   }
#'   The model will estimate a zero-inflation probability parameter \code{psi_Response} for these distributions.
#'   Example: \code{family = c(Gregarious = "binomial")}.
#' @param latent Optional character vector of latent (unmeasured) variable names.
#'   If specified, the model will account for induced correlations among observed
#'   variables that share these latent common causes.
#' @param latent_method Method for handling latent variables (default = "correlations").
#'   \itemize{
#'     \item \code{"correlations"}: MAG approach - marginalize latent variables and estimate
#'           induced correlations (\code{rho}) between observed variables that share latent parents.
#'     \item \code{"explicit"}: Model latent variables as JAGS nodes and estimate structural
#'           paths from latents to observed variables.
#'   }
#' @param standardize_latent Logical; if \code{TRUE} and \code{latent_method = "explicit"},
#'   adds standardized priors (\code{N(0,1)}) to latent variables to identify scale and location.
#'   This improves convergence and makes regression coefficients interpretable as standardized effects.
#'   Only applicable when using explicit latent variable modeling (default = TRUE).
#' @param parallel Logical; if \code{TRUE}, run MCMC chains in parallel (default = FALSE).
#'   Note: Requires \code{n.cores > 1} to take effect.
#' @param n.cores Integer; number of CPU cores to use for parallel chains (default = 1).
#'   Only used when \code{parallel = TRUE}.
#' @param cl Optional; a cluster object created by \code{parallel::makeCluster()}.
#'   If \code{NULL}, a cluster will be created and destroyed automatically.
#' @param ic_recompile Logical; if \code{TRUE} and \code{parallel = TRUE}, recompile the model
#'   after parallel chains to compute DIC/WAIC (default = TRUE).
#'   This adds a small sequential overhead but enables information criteria calculation.
#' @param optimise Logical; if \code{TRUE} (default), use the optimized random effects formulation
#'   for phylogenetic models. This is significantly faster (5-10x) and more numerically stable.
#'   If \code{FALSE}, use the traditional marginal formulation (slower, but provided for comparison).
#' @param random Optional formula or list of formulas specifying global random effects
#'   applied to all equations (e.g. \code{~(1|species)}).
#' @param levels (Hierarchical Data) A named list mapping variables to their hierarchy levels.
#'   Required if \code{data} is a list of data frames (hierarchical format).
#'   Example: \code{list(individual = c("y", "x"), site = c("z"))}.
#' @param hierarchy (Hierarchical Data) Character string describing the topological ordering of levels
#'   (e.g., \code{"site > individual"}). Required for hierarchical data if not fully inferred from random effects.
#' @param link_vars (Hierarchical Data) Optional named character vector specifying variables used to link
#'   data levels (e.g. \code{c(site = "site_id")}).
#' @param fix_residual_variance Optional named vector for fixing residual variances.
#'   Useful for handling non-identified models or specific theoretical constraints.
#'   Example: \code{c(response_var = 1)}.
#' @param priors Optional named list of character strings specifying custom priors for specific parameters.
#'   Enables overriding default uninformative priors.
#'   Example: \code{list(alpha_Response = "dnorm(0, 0.001)", beta_Response_Predictor = "dnorm(1, 10)")}.
#' @param reuse_models List of previously fitted 'because' models to scan for reusable d-separation test results.
#'   If a test in the current run matches a test in a reused model (same formula), the result is copied
#'   instead of re-running JAGS. **Note**: Ensuring that the data is consistent is the user's responsibility.
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
#' @export
#' @importFrom ape vcv.phylo branching.times
#' @importFrom rjags jags.model coda.samples dic.samples jags.samples
#' @importFrom stats na.omit update formula terms setNames start var
#' @importFrom utils capture.output
#' @importFrom coda gelman.diag effectiveSize
#' @import coda
because <- function(
  equations,
  data,
  id_col = NULL,
  structure = NULL,
  tree = NULL,
  monitor = "interpretable",
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
  parallel = FALSE,
  n.cores = parallel::detectCores() - 1,
  cl = NULL,
  ic_recompile = TRUE,
  optimise = TRUE,
  random = NULL, # Global random effects applied to all equations
  levels = NULL, # Hierarchical data: variable-to-level mapping
  hierarchy = NULL, # Hierarchical data: level ordering (e.g., "site > individual")
  link_vars = NULL, # Hierarchical data: variables linking levels
  fix_residual_variance = NULL, # Optional: fix residual variance (tau_e) for specific variables
  priors = NULL, # Optional: custom priors list
  reuse_models = NULL
) {
  # --- Input Validation & Setup ---

  # Allow single formula input
  if (inherits(equations, "formula")) {
    equations <- list(equations)
  }

  if (!quiet) {
    is_list_debug <- is.list(data) && !is.data.frame(data)
  }

  # WAIC Validity Check: Conditional WAIC (optimise=FALSE) is misleading
  if (!optimise && WAIC) {
    warning(
      "WAIC is calculated using conditional likelihoods when optimise=FALSE, which produces much higher values than marginal pseudo-likelihoods and is NOT comparable to optimise=TRUE results. Disabling WAIC to avoid confusion. Please use DIC for model comparison in this mode."
    )
    WAIC <- FALSE
  }

  latent_method <- match.arg(latent_method)

  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # Handle 'structure' alias
  if (is.null(tree) && !is.null(structure)) {
    tree <- structure
  }

  # Check if both are missing (tree is NULL, structure is NULL)
  # This implies Independent Model (tree = NULL is valid for that)

  # If both are NULL, we run as independent model.

  # Tree check moved to line 165

  # Validate inputs
  # Input validation
  if (is.null(data)) {
    stop("Argument 'data' must be provided.")
  }
  if (is.null(equations)) {
    stop("Argument 'equations' must be provided.")
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

  # --- Equation Normalization (User Request) ---
  # Support psi_Species and z_Species as aliases for Species in occupancy models.
  # This ensures the linear predictor (mu_Species) is correctly linked.
  equations <- lapply(equations, function(eq) {
    resp <- as.character(eq)[2]
    if (grepl("^(psi_|z_)", resp)) {
      base_name <- sub("^(psi_|z_)", "", resp)
      # check if base name is occupancy (explicitly or via data presence)
      is_occ <- FALSE
      if (
        !is.null(family) &&
          !is.na(family[base_name]) &&
          family[base_name] == "occupancy"
      ) {
        is_occ <- TRUE
      }
      if (base_name %in% names(data)) {
        is_occ <- TRUE
      } # Potential occupancy if matrix or vector

      if (is_occ) {
        f_str <- paste(deparse(eq), collapse = " ")
        f_str <- sub(paste0("^", resp), base_name, f_str)
        return(as.formula(f_str))
      }
    }
    return(eq)
  })

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
    }
  }

  # --- Automatic Data Cleaning (Handle Character/Factor Columns) ---
  data <- preprocess_categorical_vars(
    data,
    target_vars = model_vars,
    dummy_vars = fixed_predictors, # Categorical fixed predictors need dummies
    exclude_cols = id_col,
    quiet = quiet
  )

  # --- Hierarchical Data Detection & Validation ---
  # Data is hierarchical if it's a list (not dataframe)
  is_list_data <- is.list(data) && !is.data.frame(data)
  is_hierarchical <- FALSE
  hierarchical_info <- NULL

  if (is_list_data) {
    # Get all variables from equations for auto-detection
    eq_vars <- unique(unlist(lapply(equations, all.vars)))

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
      validate_hierarchical_data(data, levels, hierarchy, link_vars)
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

      if (!quiet) {
        message("Hierarchical data structure detected: ", hierarchy)
      }
    }
  } else {
    # Single-level data - ensure levels/hierarchy not mistakenly provided
    if (!is.null(levels) || !is.null(hierarchy)) {
      warning(
        "'levels' or 'hierarchy' provided but 'data' is not hierarchical. ",
        "Ignoring these arguments."
      )
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
    global_random_terms <- parse_global_random(random, equations)
    # Combine with equation-specific terms
    random_terms <- c(random_terms, global_random_terms)

    # Deduplicate terms (avoid adding same (1|Group) twice for same response)
    # Create unique keys
    if (length(random_terms) > 0) {
      keys <- sapply(random_terms, function(x) {
        paste(x$response, x$group, sep = "|")
      })

      if (!quiet) {
        if (length(random_terms) > 0) {
          msg <- sapply(random_terms, function(x) {
            paste(x$response, x$group, sep = "|")
          })
        } else {}
      }
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
        paste(sapply(all_poly_terms, function(x) x$original), collapse = ", ")
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
      random_vars <- unique(sapply(random_terms, function(x) x$group))
      eq_vars <- unique(c(eq_vars, random_vars))
    }

    # Remove latent variables (not in data)
    if (!is.null(latent)) {
      eq_vars <- setdiff(eq_vars, latent)
    }

    # Ensure base variables for polynomials are included (JAGS computes Age^2 from Age)
    # AND derived variables are excluded (so they aren't passed as data)
    if (!is.null(all_poly_terms)) {
      base_poly_vars <- sapply(all_poly_terms, function(x) x$base_var)
      derived_poly_vars <- sapply(all_poly_terms, function(x) x$internal_name)

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

    if (optimise && is_hierarchical) {
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
    } else {
      # OLD PATH: Flatten everything (or if not optimised?)
      # Get assembled dataset with all needed variables
      data <- get_data_for_variables(
        eq_vars,
        hierarchical_info$data,
        hierarchical_info$levels,
        hierarchical_info$hierarchy,
        hierarchical_info$link_vars
      )

      # Update original_data for later use
      original_data <- data

      if (!quiet) {
        message(
          "Assembled hierarchical data (flattened): ",
          nrow(data),
          " observations, ",
          ncol(data),
          " variables"
        )
      }
    }
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

    responses <- sapply(equations, function(eq) as.character(formula(eq)[2]))

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
    if (is.null(optimise) & is.null(tree)) {
      optimise <- TRUE
    }

    if (is_hierarchical && !is.data.frame(data) && !optimise) {
      # Should not happen if assembly worked, but safety check for flattened path
      stop("Failed to assemble hierarchical data frame for random effects.")
    }

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
      if (missing(tree) || is.null(tree)) {
        # Cannot auto-format without tree to determine species order

        warning(
          "Variability 'reps' specified but no tree provided. Automatic formatting requires a tree to order species rows. Assuming data is already aggregated or user handles index mapping."
        )
      } else if (!is.null(id_col) && id_col %in% names(data)) {
        if (!quiet) {
          message(
            "Detected 'reps' variability and long-format data. Auto-formatting matrices using because_format_data()..."
          )
        }

        # Determine the tree to use (first one if list)
        use_tree <- if (inherits(tree, "multiPhylo")) tree[[1]] else tree

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
      !(is_hierarchical && optimise)
  ) {
    # Extract all variable names from fixed equations
    eq_vars <- unique(unlist(lapply(equations, all.vars)))

    # Add variables from random terms (grouping factors)
    if (length(random_terms) > 0) {
      random_vars <- unique(sapply(random_terms, function(x) x$group))
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
    if (optimise && is_hierarchical) {
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
  if (is.null(tree)) {
    # Independent or structure via "structure" argument
    if (!is.null(structure)) {
      if (is.matrix(structure)) {
        structures[["custom"]] <- structure
      } else {
        # Generic S3 structure object - use its class name as key
        class_name <- class(structure)[1]
        structures[[class_name]] <- structure
      }
    }
  } else if (is.matrix(tree)) {
    structures[["custom"]] <- tree
  } else if (is.list(tree) && !inherits(tree, "list")) {
    # It's an S3 object with list base - use class name
    class_name <- class(tree)[1]
    # Check if it's a multi-object type (contains multiple items)
    if (length(tree) > 1 && is.null(names(tree))) {
      is_multiple <- TRUE
    }
    structures[[class_name]] <- tree
  } else if (
    is.list(tree) && is.null(class(tree)) || identical(class(tree), "list")
  ) {
    # Plain list of structures
    structures <- tree
    # Check for multi-objects in the list
    for (s in structures) {
      if (is.list(s) && length(s) > 1) is_multiple <- TRUE
    }
  } else {
    # Any other S3 object (phylo, spatial_knn, etc.) - use class name
    class_name <- class(tree)[1]
    structures[[class_name]] <- tree
  }

  structure_names <- names(structures)
  if (is.null(structure_names) && length(structures) > 0) {
    structure_names <- paste0("Struct", seq_along(structures))
    names(structures) <- structure_names
  }

  # 2. Process Structures using S3 Generic
  if (length(structures) == 0) {
    # Independent Logic: Determine N from data
    # Include matrices and arrays (useful for occupancy models)
    potential_objects <- Filter(
      function(x) is.vector(x) || is.factor(x) || is.matrix(x) || is.array(x),
      data
    )
    if (length(potential_objects) > 0) {
      obj <- potential_objects[[1]]
      N <- if (is.matrix(obj) || is.array(obj)) nrow(obj) else length(obj)
    } else {
      stop(
        "Could not determine N from data. Please provide structure or vector data."
      )
    }
  } else {
    # Use S3 Generic for Processing
    for (s_name in structure_names) {
      structure_obj <- structures[[s_name]]
      # Prepare Data
      prep_res <- prepare_structure_data(
        structure_obj,
        data = data,
        optimize = optimise,
        quiet = quiet
      )

      # Merge data updates
      if (!is.null(prep_res$data_list)) {
        for (d_name in names(prep_res$data_list)) {
          data[[d_name]] <- prep_res$data_list[[d_name]]
        }
      }

      # Determine N from the processed structure
      # Look for any square matrix or 3D array in the data_list
      current_N <- NULL
      for (d_name in names(prep_res$data_list)) {
        obj <- prep_res$data_list[[d_name]]
        if (is.matrix(obj) && nrow(obj) == ncol(obj)) {
          current_N <- nrow(obj)
          break
        } else if (is.array(obj) && length(dim(obj)) == 3) {
          current_N <- dim(obj)[1]
          break
        }
      }

      # If still NULL, check for 'n' attribute (safe way)
      if (is.null(current_N)) {
        n_attr <- attr(structure_obj, "n")
        if (!is.null(n_attr) && is.numeric(n_attr)) {
          current_N <- n_attr
        }
      }

      if (!is.null(current_N)) {
        if (is.null(N)) {
          N <- current_N
        } else if (N != current_N) {
          stop(paste("Dimension mismatch in structure:", s_name))
        }
      }

      # Update structure object if validated/standardised
      if (!is.null(prep_res$structure_object)) {
        structures[[s_name]] <- prep_res$structure_object
      }
    }

    if (optimise) {
      if (is.null(N) && "N" %in% names(data)) {
        N <- if (is.list(data)) data$N[1] else data[["N"]][1]
      } # Last resort
    }
  }

  if (is.null(N) && "N" %in% names(data)) {
    N <- if (is.list(data)) data$N[1] else data[["N"]][1]
  }
  if (is.null(N)) {
    # Try to get N from data - handle both data.frame and list cases
    if (is.data.frame(data) && nrow(data) > 0) {
      N <- nrow(data)
    } else if (is.list(data) && length(data) > 0) {
      # For lists, use length of first vector element
      first_vec <- data[[1]]
      if (is.vector(first_vec) || is.factor(first_vec)) {
        N <- length(first_vec)
      } else if (is.matrix(first_vec) || is.array(first_vec)) {
        N <- nrow(first_vec)
      }
    }
  }
  data$N <- N
  # Only add 'zeros' vector if using ZIP or ZINB (Poisson trick)
  # OR if we have structures (which often use dmnorm(zeros, ...))
  if (!is.null(N)) {
    # Check if 'zeros' already exists (use exact name match)
    has_zeros <- "zeros" %in% names(data)
    if (!has_zeros) {
      needs_zeros <- any(family %in% c("zip", "zinb", "occupancy")) ||
        length(structures) > 0

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
      response_vars <- unique(sapply(equations, function(eq) all.vars(eq[[2]])))
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

  response_vars <- unique(sapply(equations, function(eq) all.vars(eq[[2]])))
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

    # Priority 0: Occupancy Model -> Force "reps" (matrix/long format required)
    if (
      !is.null(family) &&
        var %in% names(family) &&
        family[[var]] == "occupancy"
    ) {
      auto_variability[[var]] <- "reps"
      if (!quiet) {
        message(sprintf(
          "Auto-detected: '%s' is an occupancy variable -> using 'reps' mode for detection history",
          var
        ))
      }
      next # Skip further checks
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

      # [NEW 2025-12-21] Remove occupancy variables from potential_latents for d-sep tests
      # Generic logic above adds vars not in data to 'potential_latents'.
      # Since 'SpeciesA' is renamed to 'SpeciesA_obs' in data (for occupancy models),
      # it gets added to potential_latents. We must remove it so d-sep treats it as observed.
      if (!is.null(family)) {
        occ_vars <- names(family)[family == "occupancy"]
        if (length(occ_vars) > 0) {
          potential_latents <- setdiff(potential_latents, occ_vars)
          # Also remove p_ variables if present
          p_vars <- paste0("p_", occ_vars)
          potential_latents <- setdiff(potential_latents, p_vars)
          # Also remove psi_ variables
          psi_vars <- paste0("psi_", occ_vars)
          potential_latents <- setdiff(potential_latents, psi_vars)
        }
      }

      # Exclude polynomial internal variables (they're deterministic, not latent)
      if (!is.null(all_poly_terms)) {
        poly_internal_names <- sapply(all_poly_terms, function(x) {
          x$internal_name
        })
        potential_latents <- setdiff(potential_latents, poly_internal_names)
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
    dsep_result <- because_dsep(
      equations,
      latent = latent,
      random_terms = random_terms,
      categorical_vars = if (!is.null(attr(data, "categorical_vars"))) {
        attr(data, "categorical_vars")
      } else {
        NULL
      },
      family = family,
      poly_terms = all_poly_terms,
      quiet = quiet
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

    # [User Request] Translate tests where response is a psi_ variable
    # "when psi_Species is a dependent variable, it should be tested (probably as Species ~)"
    dsep_tests <- lapply(dsep_tests, function(eq) {
      resp <- as.character(eq)[2]
      if (grepl("^psi_", resp)) {
        base_resp <- sub("^psi_", "", resp)
        # Ensure base_resp exists and is an occupancy variable
        if (
          !is.null(family) &&
            !is.na(family[base_resp]) &&
            family[base_resp] == "occupancy"
        ) {
          # Translate response to base name
          f_str <- paste(deparse(eq), collapse = " ")
          # Use regex to replace only the response at the start
          f_str <- sub(paste0("^", resp), base_resp, f_str)

          # Re-attach test_var from original
          new_eq <- as.formula(f_str)
          attr(new_eq, "test_var") <- attr(eq, "test_var")
          return(new_eq)
        }
      }
      return(eq)
    })

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
    run_single_dsep_test <- function(i, test_eq, monitor_params) {
      if (!quiet) {
        message(paste("D-sep test eq:", deparse(test_eq)))
      }

      # Select appropriate dataset for this test
      test_data <- original_data

      if (!is.null(hierarchical_info)) {
        # Extract variables from this test equation
        test_vars <- all.vars(test_eq)

        # [FIX] Add random effect grouping variables to test_vars
        # Otherwise get_data_for_variables removes them, causing "Unknown variable N_SiteID"
        if (!is.null(random_terms) && length(random_terms) > 0) {
          random_groups <- unique(sapply(random_terms, function(x) x$group))
          test_vars <- unique(c(test_vars, random_groups))
        }

        # [FIX] Handle categorical dummy variables
        # We need their PARENT variables to fetch the data, then recreate dummies manually.
        dummies_to_create <- list()

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
          hierarchical_info$data,
          hierarchical_info$levels,
          hierarchical_info$hierarchy,
          hierarchical_info$link_vars
        )

        # [FIX] Recreate dummy variables in test_data
        for (parent_var in names(dummies_to_create)) {
          if (parent_var %in% names(test_data)) {
            dummies <- dummies_to_create[[parent_var]]
            vals <- test_data[[parent_var]]

            # Recreate dummies: sex_m = as.integer(sex == 2) etc.
            levels_map <- attr(original_data, "categorical_vars")[[
              parent_var
            ]]$levels

            for (k in 2:length(levels_map)) {
              # The dummies list generally corresponds to k=2..N
              expected_dummy <- dummies[k - 1]

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
      dsep_equations <- list(test_eq)

      if (!is.null(family) && any(family == "occupancy")) {
        # Check if this test involves occupancy variables
        # We use a loop to ensure we catch supporting equations for newly added variables (recursion)
        added <- TRUE
        while (added) {
          added <- FALSE
          current_vars <- unique(unlist(lapply(dsep_equations, all.vars)))

          for (v in current_vars) {
            # Identify variable type
            # Is v an occupancy variable? (z_X or X)
            is_occ <- !is.na(family[v]) && family[v] == "occupancy"

            # Is v a detection probability p_X?
            is_det <- grepl("^p_", v) &&
              !is.na(family[sub("^p_", "", v)]) &&
              family[sub("^p_", "", v)] == "occupancy"

            # Is v a psi probability psi_X?
            is_psi <- grepl("^psi_", v) &&
              !is.na(family[sub("^psi_", "", v)]) &&
              family[sub("^psi_", "", v)] == "occupancy"

            base_occ <- if (is_det) {
              sub("^p_", "", v)
            } else if (is_psi) {
              sub("^psi_", "", v)
            } else {
              v
            }

            if (is_occ || is_det || is_psi) {
              test_eq_resp <- as.character(test_eq)[2]

              # 1. Include detection equation (p_X ~ ...) IF it is NOT the response
              p_name <- paste0("p_", base_occ)
              if (p_name != test_eq_resp) {
                for (eq in equations) {
                  if (as.character(eq)[2] == p_name) {
                    # Check if already in dsep_equations (compare as strings)
                    eq_str <- paste(deparse(eq), collapse = " ")
                    exists <- any(sapply(dsep_equations, function(e) {
                      paste(deparse(e), collapse = " ") == eq_str
                    }))
                    if (!exists) {
                      dsep_equations <- c(dsep_equations, list(eq))
                      added <- TRUE
                    }
                  }
                }
              }

              # 2. Include base occupancy equation (X ~ ...)
              should_include_occ <- FALSE
              if (is_occ && v != test_eq_resp) {
                should_include_occ <- TRUE
              }
              if (is_det && base_occ != test_eq_resp) {
                should_include_occ <- TRUE
              }
              if (
                is_psi &&
                  base_occ != test_eq_resp &&
                  paste0("psi_", base_occ) != test_eq_resp
              ) {
                should_include_occ <- TRUE
              }

              if (should_include_occ) {
                for (eq in equations) {
                  if (as.character(eq)[2] == base_occ) {
                    eq_str <- paste(deparse(eq), collapse = " ")
                    exists <- any(sapply(dsep_equations, function(e) {
                      paste(deparse(e), collapse = " ") == eq_str
                    }))
                    if (!exists) {
                      dsep_equations <- c(dsep_equations, list(eq))
                      added <- TRUE
                    }
                  }
                }
              }
            }
          }
        }
      }

      # Now filter variability/family based on ALL variables in dsep_equations
      all_dsep_vars <- unique(unlist(lapply(dsep_equations, all.vars)))
      all_dsep_vars_clean <- unique(c(
        all_dsep_vars,
        sub("^p_", "", all_dsep_vars),
        sub("^psi_", "", all_dsep_vars),
        sub("^z_", "", all_dsep_vars)
      ))

      sub_variability <- if (!is.null(variability)) {
        variability[names(variability) %in% all_dsep_vars_clean]
      } else {
        NULL
      }
      sub_family <- if (!is.null(family)) {
        family[names(family) %in% all_dsep_vars_clean]
      } else {
        NULL
      }

      # [User Request] Auto-detect binomial for binary response
      test_resp <- as.character(test_eq)[2]
      if (!is.null(test_data[[test_resp]])) {
        vals <- unique(na.omit(test_data[[test_resp]]))
        if (length(vals) <= 2 && all(vals %in% c(0, 1))) {
          if (is.null(sub_family)) {
            sub_family <- list()
          }
          if (is.na(sub_family[test_resp])) {
            sub_family[[test_resp]] <- "binomial"
          }
        }
      }

      if (length(sub_variability) == 0) {
        sub_variability <- NULL
      }
      if (length(sub_family) == 0) {
        sub_family <- NULL
      }

      fit <- because(
        data = test_data, # Use selected dataset
        tree = tree,
        equations = dsep_equations, # Use augmented list
        monitor = monitor_params,
        n.chains = n.chains,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        DIC = FALSE, # DIC not needed for d-sep tests
        WAIC = FALSE,
        n.adapt = n.adapt,
        quiet = quiet, # Restore user quiet setting for d-sep sub-runs
        dsep = FALSE,
        variability = sub_variability,
        family = sub_family,
        fix_residual_variance = fix_residual_variance,
        latent = NULL, # D-sep tests are on observed variables only
        latent_method = latent_method,
        parallel = FALSE, # Disable nested parallelism
        n.cores = 1,
        cl = NULL,
        ic_recompile = ic_recompile,
        random = random, # Pass global random effects to d-sep tests
        levels = levels, # Pass hierarchical info
        hierarchy = hierarchy,
        link_vars = link_vars
      )

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

        # Export necessary objects to cluster
        parallel::clusterExport(
          cl,
          c(
            "original_data",
            "tree",
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
            "run_single_dsep_test",
            "dsep_tests",
            "extract_random_effects"
          ),
          envir = environment()
        )

        # Run tests in parallel
        # Only run the necessary tests
        par_results <- parallel::parLapply(
          cl,
          tests_to_run_indices,
          function(i) {
            test_eq <- dsep_tests[[i]]

            # Use "interpretable" monitoring mode.
            current_monitor <- "interpretable"

            run_single_dsep_test(i, test_eq, current_monitor) # Pass current_monitor
          }
        )

        # Merge parallel results into new_results_list
        for (j in seq_along(tests_to_run_indices)) {
          idx <- tests_to_run_indices[j]
          new_results_list[[idx]] <- par_results[[j]]
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
          run_single_dsep_test(i, test_eq, current_monitor) # Pass current_monitor
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

    for (result in results) {
      samples <- result$samples
      param_map <- result$param_map
      model_string <- result$model
      i <- result$test_index

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
      equations = equations # Store original equations for safe reuse
    )
    class(result) <- "because"
    return(result)
  }

  # Handle latent variable method
  if (!is.null(latent)) {
    latent_method <- match.arg(latent_method)

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

          # Convert formula to character for manipulation
          eq_str <- paste(deparse(eq), collapse = " ")

          # Replace categorical variable with its dummies
          pattern <- paste0("\\b", var, "\\b")
          replacement <- paste(dummies, collapse = " + ")
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

  # JAGS model code
  model_output <- because_model(
    equations = equations,
    multi.tree = is_multiple,
    variability = variability_list,
    family = family,
    vars_with_na = response_vars_with_na,
    induced_correlations = induced_cors,
    latent = latent,
    standardize_latent = standardize_latent,
    optimise = optimise,
    structure_names = structure_names,
    structures = structures,
    random_structure_names = names(random_structures),
    random_terms = random_terms,
    poly_terms = all_poly_terms,
    categorical_vars = if (!is.null(attr(data, "categorical_vars"))) {
      attr(data, "categorical_vars")
    } else {
      NULL
    },
    priors = priors,
    hierarchical_info = if (is_hierarchical && optimise) {
      hierarchical_info
    } else {
      NULL
    }
  )

  model_string <- model_output$model
  parameter_map <- model_output$parameter_map

  model_file <- tempfile(fileext = ".jg")
  writeLines(model_string, model_file)

  # If latent variables are present and this is a standard run (dsep=FALSE),
  # print the MAG structure and basis set for user verification, as requested.
  # Display MAG structure for latent variable models (non-dsep runs)
  if (!dsep && !is.null(latent) && length(latent) > 0 && !quiet) {
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
            quiet = FALSE
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
  if (
    is.character(monitor) &&
      length(monitor) == 1 &&
      monitor %in% c("interpretable", "all")
  ) {
    monitor_mode <- monitor
    monitor <- NULL # Will be auto-detected based on mode
  }

  # Default mode is "interpretable"
  if (is.null(monitor_mode) && is.null(monitor)) {
    monitor_mode <- "interpretable"
  }

  if (is.null(monitor)) {
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
        # For occupancy models, the response itself (Y) is NOT a node,
        # but z_Y is. Map accordingly.
        adj_response_vars <- sapply(response_vars_all, function(v) {
          if (
            !is.null(family) &&
              !is.null(family[[v]]) &&
              family[[v]] == "occupancy"
          ) {
            return(paste0("z_", v))
          }
          return(v)
        })
        monitor <- unique(c(monitor, adj_response_vars))
      }
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

  # Automatic inits for occupancy models (z must be 1 where observations are 1)
  occupancy_inits <- list()
  if (!is.null(family)) {
    for (var_name in names(family)) {
      if (family[[var_name]] == "occupancy") {
        # Robust initialization: z=1 everywhere
        # This is compatible with all data (y=0 or y=1) and ensures correct length
        if ("N" %in% names(data)) {
          z_init <- rep(1, data$N)
          occupancy_inits[[paste0("z_", var_name)]] <- z_init
        }
      }
    }
  }

  # Clean up data list: Remove variables not present in the model code to avoid JAGS warnings
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
    if (!quiet) {
      # message(sprintf("Removing unused variables from data: %s", paste(vars_to_remove, collapse=", ")))
    }
    for (v in vars_to_remove) {
      data[[v]] <- NULL
    }
  }

  # Run MCMC chains (parallel or sequential)

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
        occupancy_inits,
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
    if (!quiet) {
      # message("--- JAGS Model Code (Pre-Compile) ---")
      # message(paste(model_output$model, collapse = "\n"))
      # message("-------------------------------------")
    }
    # Compile model

    model <- tryCatch(
      {
        if (verbose) {
          cat("\n--- JAGS MODEL STRING ---\n", model_string, "\n")
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
            occupancy_inits,
            list(
              .RNG.name = "base::Wichmann-Hill",
              .RNG.seed = 12345 + i
            )
          )
        })

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

  # Check if we can store species/unit identifiers for easy reference
  if (!is.null(tree)) {
    # If tree used, data is sorted by tip labels
    if (inherits(tree, "multiPhylo")) {
      result$species_order <- tree[[1]]$tip.label
    } else {
      result$species_order <- tree$tip.label
    }
  } else if (!is.null(id_col) && is.data.frame(original_data)) {
    # If no tree but ID col provided
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
        occupancy_inits,
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
      if (n.iter > n.burnin) {
        result$DIC <- rjags::dic.samples(model, n.iter = n.iter - n.burnin)
      } else {
        result$DIC <- NULL
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

  return(result)
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
  quiet = FALSE
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
          quiet = TRUE # Suppress output for recursive calls to avoid noise
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

  char_cols <- sapply(data[check_cols], function(x) {
    is.character(x) || is.factor(x)
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
      # Convert to factor first to get levels
      f_vals <- factor(data[[col]])
      levels <- levels(f_vals)

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
        categorical_vars[[col]] <- list(
          levels = levels,
          reference = levels[1],
          dummies = paste0(col, "_", levels[-1])
        )

        # Convert to integer codes for JAGS
        data[[col]] <- as.integer(f_vals)

        # Generate Dummy Variables explicitly ONLY if requested
        # This is required so JAGS can find 'sex_m' etc.
        # We only do this for variables used as fixed predictors to save memory.
        if (is.null(dummy_vars) || col %in% dummy_vars) {
          if (!quiet && length(levels) > 500) {
            message(sprintf(
              "Generating %d dummy variables for '%s'... this may take a moment.",
              length(levels) - 1,
              col
            ))
          }
          for (k in 2:length(levels)) {
            lev_name <- levels[k]
            dummy_col_name <- paste0(col, "_", lev_name)
            # Create binary column
            data[[dummy_col_name]] <- as.numeric(data[[col]] == k)
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
