#' Generate a JAGS model string for Bayesian SEM (Because)
#'
#' This function builds the model code to be passed to JAGS based on a set of structural equations.
#' It supports custom covariance structures (spatial, phylogenetic, etc.).
#' Missing values are handled both in the response and predictor variables treating all of them as stochastic nodes.
#'
#' @param equations A list of model formulas.
#' @param is_multi_structure Logical; if \code{TRUE}, handles 3D Precision arrays (e.g. multi-object sampling).
#' @param variability Optional character vector or named character vector of variable names that have measurement error or within-species variability.
#'   If named, the names should be the variable names and the values should be the type of variability: "se" (for mean and standard error) or "reps" (for repeated measures).
#'   If unnamed, it defaults to "se" for all specified variables.
#'   \itemize{
#'     \item "se": Expects \code{Var_mean} and \code{Var_se} in the data. The model fixes observation error: \code{Var_mean ~ dnorm(Var, 1/Var_se^2)}.
#'     \item "reps": Expects \code{Var_obs} (matrix) and \code{N_reps_Var} (vector) in the data. The model estimates observation error: \code{Var_obs[i,j] ~ dnorm(Var[i], Var_tau)}.
#'   }
#' @param family Optional named character vector specifying the family/distribution for response variables.
#'   Default is "gaussian" for all variables. Supported values: "gaussian", "binomial", "multinomial".
#'   For "binomial" variables, the model uses a logit link and a Bernoulli likelihood, with phylogenetic correlation modeled on the latent scale.
#' @param vars_with_na Optional character vector of response variable names that have missing data.
#'   These variables will use element-wise likelihoods instead of multivariate normal.
#' @param induced_correlations Optional list of variable pairs with induced correlations
#'   from latent variables. Each element should be a character vector of length 2 specifying
#'   the pair of variables that share a latent common cause.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{model}: A character string containing the JAGS model code.
#'   \item \code{parameter_map}: A data frame mapping response variables to their predictors and parameter names.
#' }
#'
#' @details
#' The generated model includes:
#' \itemize{
#'   \item Linear predictors and multivariate normal likelihoods for each response variable.
#'   \item Priors for intercepts (\code{alpha}), slopes (\code{beta}), and residual precisions (\code{tau}).
#'   \item Custom covariance modeled via provided structural objects (e.g. VCV matrices).
#'   \item (Optional) Observation models for variables with measurement error:
#'     \itemize{
#'       \item Type "se": \code{Var_mean ~ dnorm(Var, 1/Var_se^2)}
#'       \item Type "reps": \code{Var_obs[i,j] ~ dnorm(Var[i], Var_tau)}
#'     }
#'   \item (Optional) Generalized linear mixed models for non-Gaussian responses (e.g., binomial).
#'   \item (Optional) Element-wise likelihoods for response variables with missing data.
#' }
#'
#' @examples
#' eqs <- list(BR ~ BM, S ~ BR, G ~ BR, L ~ BR)
#' cat(because_model(eqs, is_multi_structure = TRUE)$model)
#'
#' @param standardize_latent Logical (default TRUE). If TRUE, standardizes latent variables to unit variance.
#' @param structures A named list of structural objects (e.g. matrices, trees) to include as correlations.
#' @param is_multi_structure Logical (Internal). If TRUE, handles 3D Precision arrays.
#' @param latent Optional character vector of latent variable names.
#' @param poly_terms (Internal) List of polynomial terms for model generation.
#' @param fix_residual_variance Optional numeric value or named vector to fix residual variance.
#' @param latent_method Method for handling latent variables ("correlations" or "explicit").
#' @param random_structure_names Optional character vector of structural names applied to all variables.
#' @param random_terms Optional list of random effects (response, group).
#' @param categorical_vars Optional character vector of categorical variable names.
#' @param priors Optional named list of custom priors.
#' @param hierarchical_info Optional list containing data hierarchy (levels, link_vars).
#' @param engine Bayesian engine to use ("jags" or "nimble").
#' @export
#' @importFrom stats formula terms setNames sd
#' @importFrom utils combn
#'
because_model <- function(
  equations,
  is_multi_structure = FALSE,
  latent_method = "correlations",
  structures = list(),
  random_structure_names = NULL,
  random_terms = list(),
  vars_with_na = NULL,
  induced_correlations = NULL,
  variability = NULL,
  family = NULL,
  standardize_latent = TRUE,
  poly_terms = NULL,
  latent = NULL,
  categorical_vars = NULL,
  fix_latent = "loading",
  fix_residual_variance = NULL,
  priors = NULL,
  hierarchical_info = NULL,
  engine = "jags",
  quiet = FALSE
) {
  # Standardize structures
  if (is.null(structures)) {
    structures <- list()
  }
  structure_names <- names(structures)
  
  # --- INTERNAL HELPERS (Toolbox) ---
  
  # Forward-compatibility: multi.tree used to be a standalone flag
  if (!is_multi_structure && "multi.tree" %in% names(match.call())) {
     is_multi_structure <- TRUE
  }
  
  # Global JAGS node tracker to prevent redefinition errors
  # Used to de-duplicate priors/derivations across different logical blocks
  declared_nodes <- character(0)
  
  # Helper: Add lines to JAGS model while preventing node redefinition
  safe_add_lines <- function(current_lines, new_lines) {
    if (length(new_lines) == 0) return(current_lines)
    
    clean_lines <- c()
    for (line in new_lines) {
      trimmed <- trimws(line)
      if (trimmed == "" || startsWith(trimmed, "#")) {
        clean_lines <- c(clean_lines, line)
        next
      }
      
      # Identify the target node (left side of ~ or <-)
      # 1. Clean up line for splitting
      # 2. Split by JAGS operators ~ or <-
      # 3. Take the left part, trim whitespace, and strip indexing [...]
      if (grepl("[~]|<-", trimmed)) {
        # Split by the first occurrence of ~ or <-
        parts <- strsplit(trimmed, "(\\~|<-)")[[1]]
        if (length(parts) > 1) {
          # The node is the first part (strip indexing brackets)
          target <- trimws(parts[1])
          target <- sub("\\[.*\\]", "", target)
          
          if (target %in% declared_nodes) {
            # SILENT FILTER: Skip redundant declarations to ensure a clean model string.
            # This allows structural hooks and prior loops to share unified names.
            next
          }
          declared_nodes <<- c(declared_nodes, target)
        }
      }
      clean_lines <- c(clean_lines, line)
    }
    return(c(current_lines, clean_lines))
  }
  
  # Helper: returns b if a is NULL or if a is a list element that doesn't exist
  `%||%` <- function(a, b) {
    tryCatch(if (!is.null(a)) a else b, error = function(e) b)
  }

  # Helper: Get prior for a parameter (custom override or default)
  get_prior <- function(param_name, type = "beta", default = NULL) {
    # 1. Check for manual override in the 'priors' argument
    if (!is.null(priors) && param_name %in% names(priors)) {
      return(paste0(param_name, " ~ ", priors[[param_name]]))
    }
    
    # 2. Check for global type override (e.g. prior = list(beta = "dnorm(...)"))
    if (!is.null(priors) && type %in% names(priors)) {
        return(paste0(param_name, " ~ ", priors[[type]]))
    }
    
    # 3. Robust Defaults (Weakly Informative SD=10 / Prec=0.01)
    # This fixes Rhat=1.8 issues while staying uninformative for human-scale data.
    default_prior <- default %||% "dnorm(0, 0.01)"
    
    # Detect if the prior is actually a fixed constant (e.g. "1.0")
    # Fixed constants use the <- operator in JAGS
    if (grepl("^[0-9.]+$", trimws(default_prior))) {
        return(paste0(param_name, " <- ", default_prior))
    }
    
    return(paste0(param_name, " ~ ", default_prior))
  }

  # Helper: Get precision prior by family
  get_precision_prior <- function(param_name, var_name) {
    if (!is.null(priors) && param_name %in% names(priors)) {
      return(paste0(param_name, " ~ ", priors[[param_name]]))
    }
    var_family <- if (!is.null(family) && var_name %in% names(family)) {
      family[[var_name]]
    } else {
      "gaussian"
    }
    fam_obj <- get_family(var_family)
    return(jags_family_precision_prior(fam_obj, param_name))
  }

  # Helper: Get hierarchical level for a variable
  get_var_level <- function(var, h_info) {
    if (is.null(h_info)) {
      return(NULL)
    }
    # Case-insensitive level lookup
    lvl_names <- names(h_info$levels)
    for (lvl_name in lvl_names) {
      # IDENTITY CHECK: Is the 'var' the level name itself?
      if (tolower(var) == tolower(lvl_name)) {
        return(lvl_name)
      }
      # MEMBERSHIP CHECK: Is the 'var' a variable within this level?
      if (any(tolower(var) == tolower(h_info$levels[[lvl_name]]))) {
        return(lvl_name)
      }
    }

    # If not found directly, check if it's a dummy variable generated from a parent categorical variable
    if (!is.null(categorical_vars)) {
      for (cat_name in names(categorical_vars)) {
        if (
          var %in%
            categorical_vars[[cat_name]]$dummies ||
            var == paste0(cat_name, "_dummy")
        ) {
          # Found the parent variable. Recursively lookup tracking the parent's level.
          return(get_var_level(cat_name, h_info))
        }
      }
    }

    return(NULL) # Variable not found in any level
  }

  # Helper: Get hierarchical level name for a random effect grouping variable
  get_random_level <- function(response, group_var, h_info) {
    if (is.null(h_info)) {
      return(NULL)
    }
    return(get_var_level(group_var, h_info))
  }

  # Helper: Get the finest level that is a descendant of all given variables' levels
  get_finest_level_of_vars <- function(vars, h_info) {
    if (is.null(h_info) || is.null(h_info$hierarchy)) {
      return(NULL)
    }

    v_lvls <- unique(na.omit(sapply(vars, function(v) {
      get_var_level(v, h_info)
    })))
    if (length(v_lvls) == 0) {
      return(NULL)
    }
    if (length(v_lvls) == 1) {
      return(v_lvls)
    }

    # Find all levels that are descendants of ALL v_lvls
    paths <- strsplit(h_info$hierarchy, "\\s*;\\s*")[[1]]
    path_list <- lapply(paths, function(p) {
      trimws(strsplit(p, "\\s*>\\s*")[[1]])
    })
    all_lvls <- unique(unlist(path_list))

    candidates <- c()
    for (lvl in all_lvls) {
      match_all <- TRUE
      for (v_lvl in v_lvls) {
        if (!is_valid_structure_mapping(v_lvl, lvl, h_info, allow_identity = TRUE)) {
          match_all <- FALSE
          break
        }
      }
      if (match_all) {
        candidates <- c(candidates, lvl)
      }
    }

    if (length(candidates) == 0) {
      return(NULL)
    }

    # Return the finest of the descendants (the intersection)
    finest <- candidates[1]
    if (length(candidates) > 1) {
      for (i in 2:length(candidates)) {
        if (is_valid_structure_mapping(finest, candidates[i], h_info, allow_identity = TRUE)) {
          finest <- candidates[i]
        }
      }
    }
    return(finest)
  }

  # Helper: Get level for a structure
  get_struct_lvl <- function(s_name, h_info) {
    if (is.null(h_info) || is.null(h_info$structure_levels)) {
      return(NULL)
    }
    return(h_info$structure_levels[[s_name]])
  }

  # Helper: Get loop bound (N or N_<Level>) for a response variable
  get_loop_bound <- function(response, h_info) {
    if (is.null(h_info)) {
      return("N")
    }
    lvl <- get_var_level(response, h_info)
    if (is.null(lvl)) {
      return("N")
    } # Fallback to N
    return(paste0("N_", lvl))
  }

  # Helper: Get zeros vector name (zeros or zeros_<Level>)
  get_zeros_name <- function(response, h_info) {
    if (is.null(h_info)) {
      return("zeros")
    }
    lvl <- get_var_level(response, h_info)
    if (is.null(lvl)) {
      return("zeros")
    }
    return(paste0("zeros_", lvl))
  }

  # Helper: Get index expression for accessing a predictor from a coarser level
  get_pred_index <- function(pred, response_level, h_info) {
    if (is.null(h_info)) {
      return(paste0(pred, "[i]"))
    }

    pred_lvl <- get_var_level(pred, h_info)
    if (is.null(pred_lvl)) {
      return(paste0(pred, "[i]"))
    } # Unknown, use default

    if (is.null(response_level) || pred_lvl == response_level) {
      return(paste0(pred, "[i]"))
    }

    # Predictor is from a coarser level
    # We need the index vector: <PredLevel>_idx_<RespLevel>
    idx_name <- paste0(pred_lvl, "_idx_", response_level)
    return(paste0(pred, "[", idx_name, "[i]]"))
  }

  # Helper: Validate if a structure's level can map to a response's level
  is_valid_structure_mapping <- function(s_lvl, r_lvl, h_info, allow_identity = FALSE) {
    # [HIERARCHY FIX] Block random effects at the SAME level as response (Identifiability)
    # UNLESS they are structured covariance components (Phylo/Spatial) where partitioning is valid.
    if (is.null(s_lvl) || is.null(r_lvl)) {
      return(TRUE)
    }
    if (s_lvl == r_lvl && !allow_identity) {
      return(FALSE)
    }
    if (is.null(h_info) || is.null(h_info$hierarchy)) {
      return(TRUE)
    }

    paths <- strsplit(h_info$hierarchy, "\\s*;\\s*")[[1]]
    for (path in paths) {
      levels <- trimws(strsplit(path, "\\s*>\\s*")[[1]])
      s_idx <- match(s_lvl, levels)
      r_idx <- match(r_lvl, levels)
      # Structure must be at or above response in the hierarchy
      if (!is.na(s_idx) && !is.na(r_idx) && s_idx <= r_idx) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  # Helper: Get index for mapping structure levels to response levels
  get_struct_index <- function(s_name, response, h_info) {
    s_lvl <- get_struct_lvl(s_name, h_info)
    r_lvl <- get_var_level(response, h_info)
    
    # [ROBUSTNESS FIX] Ensure we never return an empty index
    if (is.null(s_lvl) || is.null(r_lvl) || s_lvl == r_lvl) {
      return("i")
    }
    
    # If levels are valid but don't match, use the bridge index.
    # We verify the level names to prevent empty brackets like _idx_obs[i]
    if (nchar(s_lvl) == 0 || nchar(r_lvl) == 0) return("i")
    
    return(paste0(s_lvl, "_idx_", r_lvl, "[i]"))
  }

  # Helper: Check if random effect grouping level is compatible with response level
  is_valid_random_level <- function(response, group_var, h_info) {
    if (is.null(h_info)) {
      return(TRUE)
    }

    resp_lvl <- get_var_level(response, h_info)
    grp_lvl <- get_var_level(group_var, h_info)

    if (is.null(resp_lvl) || is.null(grp_lvl)) {
      # If hierarchy exists, but levels are unknown, assume invalid (blocked) for safety
      if (!is.null(h_info$hierarchy)) { return(FALSE) }
      return(TRUE) 
    }

    # Use the robust checker
    return(is_valid_structure_mapping(grp_lvl, resp_lvl, h_info))
  }

  # Helper: Resolve group indexing for hierarchical random effects
  get_group_idx_string <- function(response, r_name, hierarchical_info) {
    # Default behavior: group_var[i]
    group_idx <- paste0("group_", r_name, "[i]")

    if (!is.null(hierarchical_info)) {
      resp_bound <- get_loop_bound(response, hierarchical_info)
      grp_bound <- get_loop_bound(r_name, hierarchical_info)

      # Strip N_ prefix
      resp_lvl <- sub("^N_", "", resp_bound)
      grp_lvl <- sub("^N_", "", grp_bound)

      # Resolve "N" to finest level name
      finest_lvl <- trimws(strsplit(
        hierarchical_info$hierarchy,
        ">"
      )[[1]])
      finest_lvl <- finest_lvl[length(finest_lvl)]

      if (resp_lvl == "N") {
        resp_lvl <- finest_lvl
      }
      if (grp_lvl == "N") {
        grp_lvl <- finest_lvl
      }

      # If levels differ, use nested transitive index (bridge)
      if (resp_lvl != grp_lvl) {
        group_idx <- paste0(
          "group_",
          r_name,
          "[",
          grp_lvl,
          "_idx_",
          resp_lvl,
          "[i]]"
        )
      }
    }
    return(group_idx)
  }

  # Helper: Get finest level for main loop N
  main_loop_N <- "N"
  if (!is.null(hierarchical_info)) {
    h_levels <- trimws(strsplit(hierarchical_info$hierarchy, ">")[[1]])
    finest_level <- h_levels[length(h_levels)]
    main_loop_N <- paste0("N_", finest_level)
  }

  # Helper: Get N string (for hierarchical-aware array dimensions)
  N_str <- function() main_loop_N

  # Helper: Check if a structure is multi-object
  is_struct_multi <- function(s_name) {
    # 1. Consult hierarchy if available
    if (!is.null(hierarchical_info$structure_multi) && s_name %in% names(hierarchical_info$structure_multi)) {
      status <- hierarchical_info$structure_multi[[s_name]]
      if (!is.null(status)) return(as.logical(status))
    }

    # 2. Case-specific fallback for common names
    # 2. Case-specific fallback (extensions provide these via hierarchical_info or S3)
    return(FALSE)
  }

  # [HARD GATE] Uncompromising Deduplication: Purge random effects handled by structures
  if (!is.null(random_structure_names) && length(structure_names) > 0) {
    # 1. Block by Name (Immediate collision)
    random_structure_names <- setdiff(random_structure_names, structure_names)
    
    # 2. Block by Hierarchical Level (Deep collision)
    if (!is.null(hierarchical_info)) {
      struct_lvls <- tolower(unique(na.omit(sapply(structure_names, function(s) get_struct_lvl(s, hierarchical_info)))))
      if (length(struct_lvls) > 0) {
        keep_r <- vapply(random_structure_names, function(r) {
          # [NUGGET FIX] Allow random intercepts to co-exist with structured covariance 
          # at the same level (e.g. spatial Site + random Site). 
          # We only filter if the NAMES are identical (handled by setdiff above).
          return(TRUE) 
        }, logical(1))
        random_structure_names <- random_structure_names[keep_r]
      }
    }
  }

  has_structure <- !is.null(structure_names) && length(structure_names) > 0
  has_random <- !is.null(random_structure_names) &&
    length(random_structure_names) > 0
  independent <- !has_structure && !has_random

  # Flag for multi-structure (e.g., 3D precision arrays)
  # This applies to ANY structure type when is_multi_structure is TRUE

  beta_counter <- list()
  response_counter <- list()

  # Parse equations and extract deterministic terms
  deterministic_terms <- extract_deterministic_terms(equations) # Requires R/deterministic_nodes.R logic

  eq_list <- lapply(equations, function(eq) {
    response <- as.character(formula(eq))[2]
    predictors <- attr(terms(formula(eq)), "term.labels")

    # Sanitize predictors (e.g. A:B -> A_x_B)
    predictors <- sapply(predictors, sanitize_term_name)
    names(predictors) <- NULL # Remove names from sapply output

    list(response = response, predictors = predictors)
  })

  # Track all variables (for imputation)
  all_vars <- unique(unlist(lapply(eq_list, function(eq) {
    c(eq$response, eq$predictors)
  })))

  # Identify variables involved in induced correlations
  correlated_vars <- if (!is.null(induced_correlations)) {
    unique(unlist(induced_correlations))
  } else {
    character(0)
  }

  # Handle variability argument
  variability_list <- list()
  if (!is.null(variability)) {
    if (is.null(names(variability))) {
      # Default to "se" if unnamed
      variability_list <- setNames(rep("se", length(variability)), variability)
    } else {
      variability_list <- variability
    }
  }

  # Handle family argument
  dist_list <- list()
  if (!is.null(family)) {
    dist_list <- as.list(family)
  }
  param_map <- list()
  vars_error_terms <- list() # Track variables and their error terms (MAG, etc.)

  # --- Setup Common JAGS Structures ---
  # Check if we need zero_vec (for induced_correlations, multinomial, or structures)
  need_zero_vec <- !is.null(induced_correlations) || !is.null(structures)
  if (!need_zero_vec && !is.null(family)) {
    if (any(family == "multinomial") || any(family == "ordinal")) {
      need_zero_vec <- TRUE
    }
  }

  model_lines <- c(
    "model {",
    "  # Common structures and priors"
  )

  # NOTE: zero_vec and ID2 are now provided in the data list by the because() caller for reliability.
  # We do NOT define them here to avoid multiple definition errors in JAGS.

  # --- Handle Induced Correlations (MAG) Pair Processing ---
  if (!is.null(induced_correlations)) {
    model_lines <- c(
      model_lines,
      "  # Induced Correlations (Latent Variables) - Pair Processing"
    )

    # 1. Process pairs to generate correlated error terms
    for (pair in induced_correlations) {
      var1 <- pair[1]
      var2 <- pair[2]

      # Initialize lists if needed
      if (is.null(vars_error_terms[[var1]])) {
        vars_error_terms[[var1]] <- c()
      }
      if (is.null(vars_error_terms[[var2]])) {
        vars_error_terms[[var2]] <- c()
      }

      # Define pair-specific names
      res_err <- paste0("err_res_", var1, "_", var2)
      tau_res_matrix <- paste0("TAU_res_", var1, "_", var2)
      cov_matrix <- paste0("cov_", var1, "_", var2)

      model_lines <- c(
        model_lines,
        paste0(
          "  # Correlated residuals between ",
          var1,
          " and ",
          var2,
          " (Wishart Prior)"
        ),
        paste0("  ", tau_res_matrix, "[1:2, 1:2] ~ dwish(ID2[1:2, 1:2], 3)"),
        paste0(
          "  ",
          cov_matrix,
          "[1:2, 1:2] <- inverse(",
          tau_res_matrix,
          "[1:2, 1:2])"
        ),
        paste0(
          "  sigma_res_",
          var1,
          "_",
          var2,
          " <- sqrt(",
          cov_matrix,
          "[1, 1])"
        ),
        paste0(
          "  sigma_res_",
          var2,
          "_",
          var1,
          " <- sqrt(",
          cov_matrix,
          "[2, 2])"
        ),
        paste0(
          "  rho_",
          var1,
          "_",
          var2,
          " <- ",
          cov_matrix,
          "[1, 2] / (sigma_res_",
          var1,
          "_",
          var2,
          " * sigma_res_",
          var2,
          "_",
          var1,
          ")"
        )
      )

      # Generate correlated error terms from MVN
      loop_bound <- get_loop_bound(var1, hierarchical_info)
      # [SYNTAX GUARD] Ensure loop_bound is NEVER empty
      if (is.null(loop_bound) || nchar(trimws(loop_bound)) == 0) loop_bound <- "N"
      
      model_lines <- c(
        model_lines,
        paste0("  for (i in 1:", loop_bound, ") {"),
        paste0(
          "    ",
          res_err,
          "[i, 1:2] ~ dmnorm(zero_vec[1:2], ",
          tau_res_matrix,
          "[1:2, 1:2])"
        ),
        paste0("  }")
      )

      # Record error terms for each variable
      vars_error_terms[[var1]] <- c(
        vars_error_terms[[var1]],
        paste0(res_err, "[i, 1]")
      )
      vars_error_terms[[var2]] <- c(
        vars_error_terms[[var2]],
        paste0(res_err, "[i, 2]")
      )
    }
  }

  model_lines <- safe_add_lines(model_lines, "  # Structural equations")

  # --- Generic Structure Setup ---
  # [UNIFICATION] Moved to equation-specific blocks to prevent parameter duplication.
  # Global setup code for structures (e.g., priors for OU parameters) is still needed.
  if (!is.null(structures)) {
    for (s_name in names(structures)) {
      def <- jags_structure_definition(
        structures[[s_name]],
        variable_name = "err",
        s_name = s_name,
        optimize = TRUE
      )
      if (!is.null(def$setup_code)) {
        # Only harvest setup_code (constants/priors) once here
        model_lines <- safe_add_lines(model_lines, def$setup_code)
      }
    }
  }

  # --- Handle Exogenous Latent Variables ---
  # These are variables with variability (latent) but NOT response variables (no equation)
  # JAGS needs a prior for them (e.g. X[i] ~ dnorm(0, 1.0E-06))
  # Otherwise we get "Unknown variable" error
  if (!is.null(variability)) {
    model_responses <- names(response_counter) # Will be empty here, need to infer from equations
    # Re-extract all responses from equations strictly
    all_responses <- unique(sapply(equations, function(eq) {
      as.character(formula(eq))[2]
    }))

    vars_with_variability <- names(variability_list)
    exogenous_vars <- setdiff(vars_with_variability, all_responses)

    # Legacy Multi-object Setup (Fallback if structures not provided)
    if (is_multi_structure && is.null(structures)) {
      model_lines <- c(
        model_lines,
        "  # Multi-object sampling",
        "  K ~ dcat(p_obj[])",
        "  for (k in 1:Nobj) {",
        "    p_obj[k] <- 1/Nobj",
        "  }"
      )
    }

    if (length(exogenous_vars) > 0) {
      model_lines <- c(
        model_lines,
        "  # Priors for exogenous latent variables (variable with error but no parent)",
        "  for (i in 1:",
        if (!is.null(hierarchical_info)) hierarchical_info$counts[[1]] else "N",
        ") {"
      )
      for (ex_var in exogenous_vars) {
        # Uninformative prior for the latent true value
        model_lines <- c(
          model_lines,
          paste0("    ", ex_var, "[i] ~ dnorm(0, 1.0E-06)")
        )
      }
      model_lines <- c(model_lines, "  }")
    }
  }

  # main_loop_N already computed above
  # Removed global loop opening to support hierarchical loops per equation
  # model_lines <- c(model_lines, paste0("  for (i in 1:", main_loop_N, ") {"))

  # Add polynomial transformations as deterministic nodes
  # Generated with specific loops for each term's level
  if (!is.null(poly_terms)) {
    if (length(poly_terms) > 0) {
      model_lines <- c(
        model_lines,
        "    # Deterministic polynomial transformations"
      )

      for (pt in poly_terms) {
        # Determine loop bound for this specific term
        layout_N <- get_loop_bound(pt$base_var, hierarchical_info)

        # Generate wrapped code
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:", layout_N, ") {"),
          sprintf(
            "    %s[i] <- %s[i]^%d",
            pt$internal_name,
            pt$base_var,
            pt$power
          ),
          "  }"
        )
      }
    }
  }

  # Add generic deterministic nodes (e.g., interactions, logic)
  # This generalizes the polynomial logic
  if (!is.null(deterministic_terms) && length(deterministic_terms) > 0) {
    if (length(deterministic_terms) > 0) {
      model_lines <- c(
        model_lines,
        "    # Deterministic nodes (Interactions / Logic)"
      )

      for (dt_name in names(deterministic_terms)) {
        dt <- deterministic_terms[[dt_name]]

        # Determine loop bound. For interactions, use finest level of involved vars.
        # Fallback to main_loop_N if complex.

        # Use simple regex to find all variables in the dt$expression (like A[i] * B[i])
        # We need the ones BEFORE the [i]
        vars_in_expr_raw <- unique(gsub(
          "\\[i\\]",
          "",
          unlist(regmatches(
            dt$expression,
            gregexpr("\\b\\w+\\[i\\]", dt$expression)
          ))
        ))
        vars_in_expr <- unique(c(
          all.vars(parse(text = dt$original)),
          vars_in_expr_raw
        ))

        current_N <- main_loop_N

        # Determine the finest level involved in this expression
        finest_level_of_expr <- get_finest_level_of_vars(
          vars_in_expr,
          hierarchical_info
        )

        if (!is.null(finest_level_of_expr)) {
          current_N <- paste0("N_", finest_level_of_expr)
        }

        # Rebuild the expression using hierarchical-aware indexing.
        # dt$expression was generated with plain [i] for all vars; we need to
        # replace each variable's [i] with the correct cross-level index.
        # This matters when, e.g., EOS (enviro level) appears in an interaction
        # with AgeClass (individual level) — EOS needs enviro_idx_individual[i].
        expr_hierarchical <- dt$expression
        if (!is.null(hierarchical_info) && !is.null(finest_level_of_expr)) {
          # Sort vars by length descending to avoid partial replacements
          vars_sorted <- vars_in_expr[order(
            nchar(vars_in_expr),
            decreasing = TRUE
          )]
          for (v in vars_sorted) {
            correct_idx <- get_pred_index(
              v,
              finest_level_of_expr,
              hierarchical_info
            )
            # Replace the plain v[i] with the correct hierarchical index
            expr_hierarchical <- gsub(
              paste0("\\b", v, "\\[i\\]"),
              correct_idx,
              expr_hierarchical
            )
          }
        }

        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:", current_N, ") {"),
          sprintf("    %s[i] <- %s", dt$internal_name, expr_hierarchical),
          "  }"
        )
      }
    }
  }

  # Continue with structural equations
  model_lines <- safe_add_lines(model_lines, "")

  for (j in seq_along(eq_list)) {
    eq <- eq_list[[j]]
    response <- eq$response
    predictors <- eq$predictors
    dist <- dist_list[[response]] %||% "gaussian"

    # Hierarchical Loop Wrapper
    eq_loop_N <- get_loop_bound(response, hierarchical_info)
    model_lines <- safe_add_lines(model_lines, paste0("  for (i in 1:", eq_loop_N, ") {"))

    # Check if this is a deterministic identity equation (e.g., dummy ~ I(parent == k))
    is_identity <- FALSE
    if (length(predictors) == 1 && !is.null(categorical_vars)) {
      # Is response a dummy variable?
      is_dummy <- any(sapply(categorical_vars, function(cv) {
        response %in% cv$dummies
      }))
      if (is_dummy) {
        is_identity <- TRUE
      }
    }

    if (is_identity) {
      # Deterministic assignment
      expr <- term_to_jags_expression(predictors[1])
      model_lines <- c(
        model_lines,
        paste0("    ", response, "[i] <- ", expr),
        "  }"
      )
      # No parameters to monitor or priors to add for this equation
      next
    }

    # Count and assign unique suffix for the response variable
    response_count <- response_counter[[response]] %||% 0
    response_count <- response_count + 1
    response_counter[[response]] <- response_count
    suffix <- if (response_count == 1) "" else as.character(response_count)

    alpha <- paste0("alpha_", response, suffix)
    param_map[[length(param_map) + 1]] <- list(
      response = response,
      predictor = "(Intercept)",
      parameter = alpha,
      equation_index = j,
      type = "coefficient"
    )
    linpred <- alpha
    for (pred in predictors) {
      key <- paste(response, pred, suffix, sep = "_")
      if (!key %in% names(beta_counter)) {
        # Use consistent naming: beta_Response_Predictor
        beta_name <- paste0("beta_", response, suffix, "_", pred)
        beta_counter[[key]] <- beta_name
      }
      beta_name <- beta_counter[[key]]

      # Check if predictor itself is an occupancy variable (latent state interaction)
      pred_dist <- dist_list[[pred]] %||% "gaussian"

      # Get hierarchical-aware predictor index
      resp_level <- get_var_level(response, hierarchical_info)
      pred_idx <- get_pred_index(pred, resp_level, hierarchical_info)

      if (pred_dist == "occupancy") {
        # Use latent state z_Predictor instead of Predictor
        # Replace pred in pred_idx with z_pred
        z_pred_idx <- gsub(paste0("^", pred), paste0("z_", pred), pred_idx)
        linpred <- paste0(linpred, " + ", beta_name, "*", z_pred_idx)
      } else {
        linpred <- paste0(linpred, " + ", beta_name, "*", pred_idx)
      }

      # Store in parameter map
      param_map[[length(param_map) + 1]] <- list(
        response = response,
        predictor = pred,
        parameter = beta_name,
        equation_index = j,
        type = "coefficient"
      )
    }

    if (dist == "gaussian" || dist == "occupancy" || grepl("^p_", response)) {
      mu <- paste0("mu_", response, suffix)
      model_lines <- safe_add_lines(model_lines, paste0("    ", mu, "[i] <- ", linpred))
    } else if (dist == "binomial") {
      # Binomial: logit(p) = linpred + error
      # error ~ dmnorm(0, TAU)
      # We define the error mean here as 0
      mu_err <- paste0("mu_err_", response, suffix)
      err_name <- paste0("err_", response, suffix)
      p <- paste0("p_", response, suffix)
      
      err_term <- ""
      if (err_name %in% names(vars_error_terms)) {
        err_term <- paste0(" + ", err_name, "[i]")
      }

      model_lines <- c(
        model_lines,
        paste0("    ", mu_err, "[i] <- 0"),
        paste0("    logit(", p, "[i]) <- ", linpred, err_term),
        paste0("    ", response, "[i] ~ dbern(", p, "[i])")
      )
      if (engine == "jags") {
        model_lines <- c(
          model_lines,
          paste0(
            "    log_lik_",
            response,
            suffix,
            "[i] <- logdensity.bern(",
            response,
            "[i], ",
            p,
            "[i])"
          )
        )
      }
    } else if (dist == "multinomial") {
      # Multinomial: K categories
      # L[i, 1] <- 0
      # L[i, k] <- alpha[k] + beta*X + err[i, k]

      K_var <- paste0("K_", response)
      err <- paste0("err_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    # Multinomial linear predictor for ", response),
        paste0("    L_", response, "[i, 1] <- 0"),
        paste0("    for (k in 2:", K_var, ") {")
      )

      # Linear predictor for k-th dimension
      # Note: We use array indexing [k] for parameters
      linpred_k <- paste0("alpha_", response, "[k]")

      for (pred in predictors) {
        # We force a specific name for multinomial betas to ensure array usage
        beta_name <- paste0("beta_", response, "_", pred)

        # Get hierarchical-aware predictor index
        resp_level <- get_var_level(response, hierarchical_info)
        pred_idx <- get_pred_index(pred, resp_level, hierarchical_info)
        linpred_k <- paste0(linpred_k, " + ", beta_name, "[k] * ", pred_idx)

        # Add intercepts to parameter map (only once)
        if (!paste0("alpha_", response) %in% names(beta_counter)) {
          beta_counter[[paste0("alpha_", response)]] <- TRUE
          param_map[[length(param_map) + 1]] <- list(
            response = response,
            predictor = "(Intercepts)",
            parameter = paste0("alpha_", response, "[]"),
            equation_index = j,
            type = "coefficient"
          )
        }

        # Map (only once)
        key <- paste(response, pred, suffix, sep = "_")
        if (!key %in% names(beta_counter)) {
          beta_counter[[key]] <- beta_name # Mark as used
          param_map[[length(param_map) + 1]] <- list(
            response = response,
            predictor = pred,
            parameter = paste0(beta_name, "[]"),
            equation_index = j,
            type = "coefficient"
          )
        }
      }

      # [SEMANTIC GUARD] Ensure linpred_k is not empty
      if (is.null(linpred_k) || nchar(trimws(linpred_k)) == 0) linpred_k <- "0"
      
      linpred_k <- paste0(linpred_k, " + ", err, "[i, k]")

      model_lines <- c(
        model_lines,
        paste0("      L_", response, "[i, k] <- max(-20, min(20, ", linpred_k, "))"),
        "    }"
      )

      model_lines <- c(
        model_lines,
        paste0("    # Softmax for ", response)
      )

      if (engine == "nimble") {
        # Numerically stable softmax for NIMBLE
        model_lines <- c(
          model_lines,
          paste0("    max_L_", response, "[i] <- max(L_", response, "[i, 1:", K_var, "])"),
          paste0("    for (k in 1:", K_var, ") {"),
          paste0(
            "      exp_L_",
            response,
            "[i, k] <- exp(L_",
            response,
            "[i, k] - max_L_",
            response,
            "[i])"
          ),
          "    }"
        )
      } else {
        # Standard JAGS softmax
        model_lines <- c(
          model_lines,
          paste0("    for (k in 1:", K_var, ") {"),
          paste0(
            "      exp_L_",
            response,
            "[i, k] <- exp(L_",
            response,
            "[i, k])"
          ),
          "    }"
        )
      }

      model_lines <- c(
        model_lines,
        paste0(
          "    sum_exp_L_",
          response,
          "[i] <- sum(exp_L_",
          response,
          "[i, 1:",
          K_var,
          "])"
        ),
        paste0("    for (k in 1:", K_var, ") {"),
        paste0(
          "      p_",
          response,
          "[i, k] <- (exp_L_",
          response,
          "[i, k] / sum_exp_L_",
          response,
          "[i]) * 0.9999999 + 1.0e-10"
        ),
        "    }",
        paste0(
          "    ",
          response,
          "[i] ~ dcat(p_",
          response,
          "[i, 1:",
          K_var,
          "])"
        )
      )
      if (engine == "jags") {
        model_lines <- c(
          model_lines,
          paste0(
            "    log_lik_",
            response,
            suffix,
            "[i] <- logdensity.cat(",
            response,
            "[i], p_",
            response,
            "[i, 1:",
            K_var,
            "])"
          )
        )
      }
    } else if (dist == "ordinal") {
      # Ordinal: Cumulative Logit (Proportional Odds)
      # P(Y <= k) = logit^(-1)(cutpoint[k] - eta)
      # eta = linpred + error

      K_var <- paste0("K_", response)
      err_name <- paste0("err_", response, suffix)
      eta <- paste0("eta_", response, suffix)

      err_term <- ""
      if (err_name %in% names(vars_error_terms)) {
        err_term <- paste0(" + ", err_name, "[i]")
      }

      # Linear predictor (eta)
      # Note: No intercept in eta (intercept is absorbed into cutpoints)
      linpred_no_int <- "0"
      if (!is.null(predictors) && length(predictors) > 0) {
        for (pred in predictors) {
          beta_name <- paste0("beta_", response, "_", pred)

          # Get hierarchical-aware predictor index
          resp_level <- get_var_level(response, hierarchical_info)
          pred_idx <- get_pred_index(pred, resp_level, hierarchical_info)

          linpred_no_int <- paste0(
            linpred_no_int,
            " + ",
            beta_name,
            " * ",
            pred_idx
          )

          # Map
          key <- paste(response, pred, suffix, sep = "_")
          if (!beta_name %in% names(beta_counter)) {
            beta_counter[[length(beta_counter)+1]] <- beta_name
            names(beta_counter)[length(beta_counter)] <- key
            param_map[[length(param_map) + 1]] <- list(
              response = response,
              predictor = pred,
              parameter = beta_name,
              equation_index = j,
              type = "coefficient"
            )
          }
        }
      }

      # Add cutpoints to parameter map
      param_map[[length(param_map) + 1]] <- list(
        response = response,
        predictor = "(Cutpoints)",
        parameter = paste0("cutpoint_", response, suffix),
        equation_index = j,
        type = "coefficient"
      )

      # [SEMANTIC GUARD] Ensure err is not empty
      if (is.null(err) || nchar(trimws(err)) == 0) err <- "0"

      model_lines <- c(
        model_lines,
        paste0("    # Ordinal linear predictor for ", response),
        paste0("    ", eta, "[i] <- ", linpred_no_int, " + ", err, "[i]"),

        # Cumulative probabilities
        paste0("    for (k in 1:(", K_var, "-1)) {"),
        paste0(
          "      logit(Q_",
          response,
          "[i, k]) <- cutpoint_",
          response,
          "[k] - ",
          eta,
          "[i]"
        ),
        "    }",

        # Category probabilities
        paste0("    p_", response, "[i, 1] <- Q_", response, "[i, 1]"),
        paste0("    for (k in 2:(", K_var, "-1)) {"),
        paste0(
          "      p_",
          response,
          "[i, k] <- Q_",
          response,
          "[i, k] - Q_",
          response,
          "[i, k-1]"
        ),
        "    }",
        paste0(
          "    p_",
          response,
          "[i, ",
          K_var,
          "] <- 1 - Q_",
          response,
          "[i, ",
          K_var,
          "-1]"
        ),

        # Likelihood
        paste0(
          "    ",
          response,
          "[i] ~ dcat(p_",
          response,
          "[i, 1:",
          K_var,
          "])"
        )
      )
      if (engine == "jags") {
        model_lines <- c(
          model_lines,
          paste0(
            "    log_lik_",
            response,
            suffix,
            "[i] <- logdensity.cat(",
            response,
            "[i], p_",
            response,
            "[i, 1:",
            K_var,
            "])"
          )
        )
      }
    } else if (dist == "poisson") {
      # Poisson: log(μ) = linpred + error
      # Only include err term if residuals are requested or present (handling overdispersion)
      
      mu <- paste0("mu_", response, suffix)
      err_term <- ""
      err_name <- paste0("err_", response, suffix)
      
      # [IDENTITY GUARD] Only include additive error if it's been registered
      if (err_name %in% names(vars_error_terms)) {
        err_term <- paste0(" + ", err_name, "[i]")
      }

      model_lines <- c(
        model_lines,
        paste0("    # Poisson log link for ", response),
        paste0(
          "    ",
          mu,
          "[i] <- exp(max(-20, min(10, ",
          linpred,
          err_term,
          ")))"
        ),
        paste0("    ", response, "[i] ~ dpois(", mu, "[i])")
      )
      if (engine == "jags") {
        model_lines <- c(
          model_lines,
          paste0(
            "    log_lik_",
            response,
            suffix,
            "[i] <- logdensity.pois(",
            response,
            "[i], ",
            mu,
            "[i])"
          )
        )
      }
    } else if (dist == "negbinomial") {
      # Negative Binomial: log(μ) = linpred + error
      # Y ~ NegBin(p, r) where p = r/(r+μ) and r = size parameter

      mu <- paste0("mu_", response, suffix)
      p <- paste0("p_", response, suffix)
      r <- paste0("r_", response, suffix)
      
      err_term <- ""
      err_name <- paste0("err_", response, suffix)
      
      # [IDENTITY GUARD] Only include additive error if it's been registered
      if (err_name %in% names(vars_error_terms)) {
        err_term <- paste0(" + ", err_name, "[i]")
      }

      model_lines <- c(
        model_lines,
        paste0("    # Negative Binomial log link for ", response),
        paste0(
          "    ",
          mu,
          "[i] <- exp(max(-20, min(10, ",
          linpred,
          err_term,
          ")))"
        )
      )
      if (engine == "jags") {
        model_lines <- c(
          model_lines,
          paste0("    ", p, "[i] <- ", r, " / (", r, " + ", mu, "[i])"),
          paste0("    ", response, "[i] ~ dnegbin(", p, "[i], ", r, ")")
        )
      } else {
        # NIMBLE optimization: use custom dnb for stability
        # Positional arguments: mu, r
        model_lines <- c(
          model_lines,
          paste0(
            "    ",
            response,
            "[i] ~ dnb_because(",
            mu,
            "[i], ",
            r,
            ")"
          )
        )
      }
      if (engine == "jags") {
        model_lines <- c(
          model_lines,
          paste0(
            "    log_lik_",
            response,
            suffix,
            "[i] <- logdensity.negbin(",
            response,
            "[i], ",
            p,
            "[i], ",
            r,
            ")"
          )
        )
      }
    } else if (dist == "zip") {
      # Zero-Inflated Poisson:
      # P(Y=0) = psi + (1-psi)*exp(-mu)
      # P(Y=y) = (1-psi)*dpois(y, mu) for y>0

      err <- paste0("err_", response, suffix)
      mu <- paste0("mu_", response, suffix)
      psi <- paste0("psi_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    # ZIP log link for ", response),
        paste0(
          "    ",
          mu,
          "[i] <- exp(max(-20, min(10, ",
          linpred,
          " + ",
          err,
          "[i])))"
        ),

        # Zeros trick for custom likelihood
        # log_lik terms
        paste0(
          "    lik_zero_",
          response,
          "[i] <- ",
          psi,
          " + (1-",
          psi,
          ") * exp(-",
          mu,
          "[i])"
        ),
        paste0(
          "    lik_pos_",
          response,
          "[i] <- (1-",
          psi,
          ") * exp(",
          if (engine == "jags") {
            paste0("logdensity.pois(", response, "[i], ", mu, "[i])")
          } else {
            paste0("dpois(", response, "[i], ", mu, "[i], 1)")
          },
          ")"
        )
      )

      # Select likelihood based on Y[i]
      if (engine == "jags") {
        # Zeros trick for custom likelihood
        model_lines <- c(
          model_lines,
          paste0(
            "    is_zero_",
            response,
            suffix,
            "[i] <- step(0.5 - ",
            response,
            "[i])"
          ),
          paste0(
            "    lik_",
            response,
            suffix,
            "[i] <- is_zero_",
            response,
            suffix,
            "[i] * lik_zero_",
            response,
            suffix,
            "[i] + (1 - is_zero_",
            response,
            suffix,
            "[i]) * lik_pos_",
            response,
            suffix,
            "[i]"
          ),
          paste0(
            "    log_lik_",
            response,
            suffix,
            "[i] <- log(max(1.0E-30, lik_",
            response,
            suffix,
            "[i]))"
          ),
          paste0(
            "    phi_",
            response,
            suffix,
            "[i] <- -log_lik_",
            response,
            suffix,
            "[i] + 10000"
          ),
          paste0("    zeros[i] ~ dpois(phi_", response, suffix, "[i])")
        )
      } else {
        # NIMBLE optimization: use direct dzip
        # Positional arguments: mu, psi
        model_lines <- c(
          model_lines,
          paste0(
            "    ",
            response,
            "[i] ~ dzip_because(",
            mu,
            "[i], psi_",
            response,
            suffix,
            ")"
          )
        )
      }
    } else if (dist == "zinb") {
      # Zero-Inflated Negative Binomial:
      # P(Y=0) = psi + (1-psi)*(r/(r+mu))^r
      # P(Y=y) = (1-psi)*dnegbin(y, p, r) for y>0

      err <- paste0("err_", response, suffix)
      mu <- paste0("mu_", response, suffix)
      r <- paste0("r_", response, suffix)
      p <- paste0("p_", response, suffix)
      psi <- paste0("psi_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    # ZINB log link for ", response),
        paste0(
          "    ",
          mu,
          "[i] <- exp(max(-20, min(10, ",
          linpred,
          " + ",
          err,
          "[i])))"
        ),
        paste0("    ", p, "[i] <- ", r, " / (", r, " + ", mu, "[i])")
      )

      if (engine == "jags") {
        # Zeros trick
        # Zero case: psi + (1-psi) * p^r
        model_lines <- c(
          model_lines,
          paste0(
            "    lik_zero_",
            response,
            suffix,
            "[i] <- ",
            psi,
            " + (1-",
            psi,
            ") * pow(",
            p,
            "[i], ",
            r,
            ")"
          ),

          # Positive case: (1-psi) * dnegbin(...)
          paste0(
            "    lik_pos_",
            response,
            suffix,
            "[i] <- (1-",
            psi,
            ") * exp(",
            if (engine == "jags") {
              paste0(
                "logdensity.negbin(",
                response,
                "[i], ",
                p,
                "[i], ",
                r,
                ")"
              )
            } else {
              paste0("dnegbin(", response, "[i], ", p, "[i], ", r, ", 1)")
            },
            ")"
          ),

          # Select likelihood
          paste0(
            "    is_zero_",
            response,
            suffix,
            "[i] <- step(0.5 - ",
            response,
            "[i])"
          ),
          paste0(
            "    lik_",
            response,
            suffix,
            "[i] <- is_zero_",
            response,
            suffix,
            "[i] * lik_zero_",
            response,
            suffix,
            "[i] + (1 - is_zero_",
            response,
            suffix,
            "[i]) * lik_pos_",
            response,
            suffix,
            "[i]"
          ),

          paste0(
            "    log_lik_zero_",
            response,
            suffix,
            "[i] <- log(max(1.0E-30, ",
            "lik_zero_",
            response,
            suffix,
            "[i]))"
          ),
          paste0(
            "    log_lik_pos_",
            response,
            suffix,
            "[i] <- log(max(1.0E-30, 1-",
            psi,
            ")) + ",
            if (engine == "jags") {
              paste0(
                "logdensity.negbin(",
                response,
                "[i], ",
                p,
                "[i], ",
                r,
                ")"
              )
            } else {
              paste0("dnegbin(", response, "[i], ", p, "[i], ", r, ", 1)")
            }
          ),
          paste0(
            "    log_lik_",
            response,
            suffix,
            "[i] <- is_zero_",
            response,
            suffix,
            "[i] * log_lik_zero_",
            response,
            suffix,
            "[i] + (1 - is_zero_",
            response,
            suffix,
            "[i]) * log_lik_pos_",
            response,
            suffix,
            "[i]"
          ),
          paste0(
            "    phi_",
            response,
            suffix,
            "[i] <- -log_lik_",
            response,
            suffix,
            "[i] + 10000"
          ),
          paste0(
            "    zeros_",
            response,
            "[i] ~ dpois(phi_",
            response,
            "[i])"
          )
        )
      } else {
        # NIMBLE optimization: use direct dzinb
        # Positional arguments: mu, r, psi
        model_lines <- c(
          model_lines,
          paste0(
            "    ",
            response,
            "[i] ~ dzinb_because(",
            mu,
            "[i], ",
            r,
            ", psi_",
            response,
            suffix,
            ")"
          )
        )
      }
    } else if (dist == "occupancy") {
      # No Zeros-trick or custom likelihood loop needed here
      # The occupancy likelihood is handled in the Likelihoods section.
    } else {
      stop(paste("Unknown distribution:", dist))
    }

    # Close Hierarchical Loop
    model_lines <- c(model_lines, "  }")
  }

  model_lines <- safe_add_lines(model_lines, "  # Multivariate normal likelihoods")

  # Likelihoods for responses
  for (response in unique(names(response_counter))) {
    # Skip if involved in induced correlations (handled separately for Gaussian)
    dist <- dist_list[[response]] %||% "gaussian"
    if (response %in% correlated_vars && dist == "gaussian") {
      next
    }

    dist <- dist_list[[response]] %||% "gaussian"
    # Note: occupancy and detection auxiliary equations are handled inside the loop

    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      dist <- dist_list[[response]] %||% "gaussian"
      response_var <- paste0(response, suffix)

      # Determine proper loop bound for this response variable
      loop_bound <- get_loop_bound(response, hierarchical_info)
      resp_data_idx <- paste0(response, "[i]")

      tau <- paste0("TAU", tolower(response), suffix)

      if (dist == "gaussian") {
        # Skip likelihood generation for detection probability equations (p_Response)
        # These are auxiliary equations for Occupancy models
        if (grepl("^p_", response)) {
          target_resp <- sub("^p_", "", response)
          # Check if this target is an occupancy model response
          is_occupancy_aux <- any(sapply(names(dist_list), function(n) {
            (dist_list[[n]] == "occupancy") && (n == target_resp)
          }))
          if (is_occupancy_aux) next
        }

        mu <- paste0("mu_", response, suffix)
        tau_scalar <- paste0("tau", response, suffix)

        # Check if this variable has missing data
        if (FALSE) {
          # Legacy handling for Gaussian missing data
        } else {
          if (independent) {
            # Independent Model (No random effects)
            # Y[i] ~ dnorm(mu[i], tau)
            # Note: We use tau_res for consistency across all model types.
            tau_res <- paste0("tau_res_", response, suffix)

            model_lines <- c(
              model_lines,
              paste0("  for (i in 1:", loop_bound, ") {"),
              paste0(
                "    ",
                response_var,
                "[i] ~ dnorm(",
                mu,
                "[i], ",
                tau_res,
                ")"
              )
            )
            if (engine == "jags") {
              model_lines <- c(
                model_lines,
                paste0(
                  "    log_lik_",
                  response,
                  suffix,
                  "[i] <- logdensity.norm(",
                  response_var,
                  "[i], ",
                  mu,
                  "[i], ",
                  tau_res,
                  ")"
                )
              )
            }
            model_lines <- safe_add_lines(model_lines, paste0("  }"))
          } else {
            # Optimized Random Effects Formulation (Additive)
            additive_terms <- ""
            # 1. Structures (Phylo, Spatial, etc.)
            processed_signals <- character(0)
            for (s_name in structure_names) {
              # [SINGULARITY GUARD] Ensure each structure is only added once per response
              if (s_name %in% processed_signals) next
              processed_signals <- c(processed_signals, s_name)
              
              if (
                !is_valid_structure_mapping(
                  get_struct_lvl(s_name, hierarchical_info),
                  get_var_level(response, hierarchical_info),
                  hierarchical_info,
                  allow_identity = TRUE
                )
              ) {
                next
              }

              # Structure properties
              s_lvl <- get_struct_lvl(s_name, hierarchical_info)
              s_bound <- if (is.null(s_lvl)) loop_bound else paste0("N_", s_lvl)
              s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)

              # Hierarchical bridge index (e.g. site_idx_obs[i])
              s_idx_var <- get_struct_index(s_name, response, hierarchical_info)

              # Determine if we should use partitioning (Lambda logic)
              # Condition: Gaussian variable + exactly one structure mapping to its level
              r_lvl <- get_var_level(response, hierarchical_info)
              
              # Count how many structures are valid for this specific response level
              local_struct_count <- 0
              for (sn in structure_names) {
                if (is_valid_structure_mapping(get_struct_lvl(sn, hierarchical_info), r_lvl, hierarchical_info, allow_identity = TRUE)) {
                  local_struct_count <- local_struct_count + 1
                }
              }
              
              use_partitioning <- (local_struct_count == 1) && 
                                 (get_family_object(dist_list[[response]] %||% "gaussian")$family == "gaussian") &&
                                 (s_lvl == r_lvl)

              # Call Generic to define the error vector u ~ dmnorm(...)
              s_def <- jags_structure_definition(
                structures[[s_name]],
                variable_name = response,
                s_name = s_name,
                loop_bound = s_bound,
                zeros_name = s_zeros,
                is_multi = is_struct_multi(s_name),
                i_index = s_idx_var,
                use_partitioning = use_partitioning
              )

              model_lines <- safe_add_lines(model_lines, s_def$model_lines)
              
              if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
                 additive_terms <- paste0(additive_terms, " + ", s_def$term)
              }

              # Map structural parameters (DEDUPLICATED)
              param_name <- paste0("sigma_", response, "_", s_name)
              tau_name <- paste0("tau_u_", response, "_", s_name)
              
              # Map structural parameters (DEDUPLICATED with Hierarchical Priority)
              param_name <- paste0("sigma_", response, "_", s_name)
              tau_name <- paste0("tau_u_", response, "_", s_name)
              
              # Find if this signal pair already exists
              existing_idx <- which(sapply(param_map, function(p) {
                p$response == response && p$predictor == s_name && p$type == "structure"
              }))
              
            # [PACKAGE FIX] Detect overlaps with structured covariances
            is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(u, s_name), logical(1)))
            
            if (is_unified) {
                # [UNIFICATION] Register the unified name instead of legacy
                unified_tau <- paste0("tau_u_", s_name, "_", response)
                unified_sigma <- paste0("sigma_", s_name, "_", response)
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
                
                # If partitioning was used, register lambda too
                if (isTRUE(s_def$partition_handled)) {
                    param_map[[length(param_map) + 1]] <- list(
                        response = response, 
                        predictor = s_name, 
                        parameter = paste0("lambda_", response), 
                        type = "structure"
                    )
                }
            } else {
              if (length(existing_idx) > 0) {
                # [PRIORITY FIX] If existing is legacy and current is hierarchical, overwrite
                existing_param <- param_map[[existing_idx[1]]]$parameter
                if (grepl(paste0("_", s_name, "$"), param_name) && !grepl(paste0("_", s_name, "$"), existing_param)) {
                   param_map[[existing_idx[1]]]$parameter <- tau_name
                   param_map[[existing_idx[2]]]$parameter <- param_name
                }
              } else {
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = tau_name, equation_index = NA, type = "structure")
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = param_name, equation_index = NA, type = "structure")
              }
            }
            }

            # 2. Random Effects (Grouped)
            for (r_name in random_structure_names) {
              if (!is_valid_random_level(response, r_name, hierarchical_info)) next
              
              # [UNIFICATION] Skip if handled as a unified structure
              if (r_name %in% names(structures) || any(vapply(names(structures), function(s) grepl(s, r_name), logical(1)))) {
                  is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
                  if (is_unified_r) next
              }
              
              s_suffix_r <- paste0("_", r_name)
              u_std <- paste0("u_std_", response, suffix, s_suffix_r)
              u <- paste0("u_", response, suffix, s_suffix_r)
              tau_u <- paste0("tau_u_", response, suffix, s_suffix_r)
              
              # Determine hierarchical level if available
              r_lvl <- get_random_level(response, r_name, hierarchical_info)
              # [SYNTAX GUARD] Final fallback to prevent empty [1:] indices
              n_groups <- if (!is.null(r_lvl) && nchar(r_lvl) > 0) paste0("N_", r_lvl) else if (nchar(r_name) > 0) paste0("N_", r_name) else "N"
              zeros_name <- if (!is.null(r_lvl) && nchar(r_lvl) > 0) paste0("zeros_", r_lvl) else if (nchar(r_name) > 0) paste0("zeros_", r_name) else "zeros"
              
              prec_name <- paste0("Prec_", r_name)

              model_lines <- safe_add_lines(model_lines, c(
                paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
                paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- ", u_std, "[g] / sqrt(", tau_u, ") }")
              ))
              
              group_idx <- get_group_idx_string(response, r_name, hierarchical_info)
              additive_terms <- paste0(additive_terms, " + ", u, "[", group_idx, "]")
              
              param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = tau_u, equation_index = NA, type = "structure")
              param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name), equation_index = NA, type = "structure")
            }

            # 3. Unified Observation Loop
            tau_res <- paste0("tau_res_", response, suffix)
            
            # [PARTITIONING GUARD] If the variance was already partitioned, 
            # the extension has already generated the tau_res prior.
            partition_already_handled <- FALSE
            if (!is.null(structures)) {
                if (exists("s_def") && isTRUE(s_def$partition_handled) && s_def$variable_name == response) {
                    partition_already_handled <- TRUE
                }
            }
            
            if (!partition_already_handled) {
                # [CORE FIX] Generate prior here if not handled by extension
                p_obj <- get_family_object(dist_list[[response]] %||% "gaussian")
                model_lines <- safe_add_lines(model_lines, jags_family_precision_prior(p_obj, tau_res))
            }

            model_lines <- c(
              model_lines,
              paste0("  for (i in 1:", loop_bound, ") {"),
              paste0(
                "    ",
                response_var,
                "[i] ~ dnorm(",
                mu,
                "[i]",
                additive_terms, # Sum all structural / random effects
                ", ",
                tau_res,
                ")"
              )
            )

            if (engine == "jags") {
              model_lines <- c(
                model_lines,
                paste0(
                  "    log_lik_",
                  response,
                  suffix,
                  "[i] <- logdensity.norm(",
                  response_var,
                  "[i], ",
                  mu,
                  "[i]",
                  additive_terms,
                  ", ",
                  tau_res,
                  ")"
                )
              )
            }
            model_lines <- safe_add_lines(model_lines, "  }")
          }
        }
      } else if (dist == "binomial") {
        # For binomial, the error term has the phylogenetic structure
        err <- paste0("err_", response, suffix)

        if (independent) {
          # Independent Binomial: Standard GLM (no residual error)
          model_lines <- c(
            model_lines,
            paste0("  # Independent (Standard GLM) for binomial: ", response)
          )
        } else {
          # Optimized Random Effects Formulation
          epsilon <- paste0("epsilon_", response, suffix)
          tau_res <- paste0("tau_res_", response, suffix)
          
          # Initialize accumulators
          total_u <- ""
          local_obs_code <- c()

          # 1. Structures (Phylo, Spatial, etc.)
          for (s_name in structure_names) {
            if (!is_valid_structure_mapping(get_struct_lvl(s_name, hierarchical_info), get_var_level(response, hierarchical_info), hierarchical_info, allow_identity = TRUE)) next
            
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) loop_bound else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)

            # Hierarchical bridge index (e.g. site_idx_obs[i])
            s_idx_var <- get_struct_index(s_name, response, hierarchical_info)

            s_def <- jags_structure_definition(
              structures[[s_name]],
              variable_name = response,
              s_name = s_name,
              loop_bound = s_bound,
              zeros_name = s_zeros,
              category_index = NULL,
              is_multi = is_struct_multi(s_name),
              i_index = s_idx_var
            )
            is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(u, s_name), logical(1)))
            if (!is.null(s_def)) {
              if (is_unified) {
                 unified_tau <- paste0("tau_u_", s_name, "_", response)
                 unified_sigma <- paste0("sigma_", s_name, "_", response)
                 param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                 param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
              } else {
                 # Map legacy structural parameters
                 param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                 param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
              }

              # Add Vector Prior
              model_lines <- safe_add_lines(model_lines, s_def$model_lines)
              
              # Add to observation math
              if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
                 total_u <- paste0(total_u, " + ", s_def$term)
              }
            }
          }

          # 2. Random Effects (Grouped)
          for (r_name in random_structure_names) {
            if (!is_valid_random_level(response, r_name, hierarchical_info)) next
            
            # [UNIFICATION] Skip if handled as a unified structure
            if (r_name %in% names(structures) || any(vapply(names(structures), function(s) grepl(s, r_name), logical(1)))) {
                is_unified_r <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(r_name)), logical(1)))
                if (is_unified_r) next
            }

            # [SMART DEDUPLICATION] Strictly honor manual terms (nuggets)
            relevant_term <- Filter(function(rt) rt$group == r_name && rt$response == response, random_terms)
            source_tag <- if (length(relevant_term) > 0) relevant_term[[1]]$source else "auto"

            s_suffix <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)
            
            r_lvl <- get_random_level(response, r_name, hierarchical_info)
            n_groups <- if (!is.null(r_lvl)) paste0("N_", r_lvl) else paste0("N_", r_name)
            zeros_name <- if (!is.null(r_lvl)) paste0("zeros_", r_lvl) else paste0("zeros_", r_name)
            prec_name <- paste0("Prec_", r_name)

            model_lines <- c(
              model_lines,
              paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
              paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- ", u_std, "[g] / sqrt(", tau_u, ") }")
            )
            
            group_idx <- get_group_idx_string(response, r_name, hierarchical_info)
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
            
            param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = tau_u, equation_index = NA, type = "structure")
            param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name), equation_index = NA, type = "structure")
          }

          # 3. Final Observation Loop for Binomial
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) paste0(" + ", paste(mag_terms, collapse = " + ")) else ""

          model_lines <- c(
            model_lines,
            paste0("  # Binomial error term summation: ", response),
            paste0("  for (i in 1:", loop_bound, ") {"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_res, ")"),
            paste0("    ", err, "[i] <- ", epsilon, "[i]", total_u, total_mag),
            paste0("  }")
          )
        }
      } else if (dist == "multinomial") {
            # Multinomial error terms: err[1:N, k]
            err <- paste0("err_", response, suffix)
            K_var <- paste0("K_", response)

            if (independent) {
              # Independent Multinomial: Standard GLM
              model_lines <- c(
                model_lines,
                paste0("  # Independent (Standard GLM) for multinomial: ", response),
                paste0("  for (k in 2:", K_var, ") {"),
                paste0("    for (i in 1:", loop_bound, ") {"),
                paste0("      ", err, "[i, k] <- 0"),
                paste0("    }"),
                paste0("  }")
              )
            } else {
              # Optimized Random Effects Formulation for Multinomial
              epsilon <- paste0("epsilon_", response, suffix)
              tau_res <- paste0("tau_res_", response, suffix)

              model_lines <- c(
                model_lines,
                paste0("  # Random effects for multinomial: ", response),
                paste0("  for (k in 2:", K_var, ") {")
              )

              # Initialize accumulators for category k
              total_u <- ""
              
              # 1. Custom Structures (Hooks)
              if (length(structures) > 0) {
                processed_multi_signals <- character(0)
                for (s_idx in seq_along(structures)) {
                  s_name <- names(structures)[s_idx]
                  if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
                  
                  # [SINGULARITY GUARD]
                  if (s_name %in% processed_multi_signals) next
                  processed_multi_signals <- c(processed_multi_signals, s_name)
                  
                  s_obj <- structures[[s_idx]]

                  # [LINEAGE GUARD] Only process valid hierarchical branches
                  if (!is_valid_structure_mapping(get_struct_lvl(s_name, hierarchical_info), get_var_level(response, hierarchical_info), hierarchical_info, allow_identity = TRUE)) next

                  s_lvl <- get_struct_lvl(s_name, hierarchical_info)
                  s_bound <- if (is.null(s_lvl)) get_loop_bound(response, hierarchical_info) else paste0("N_", s_lvl)
                  s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)

                  # Hierarchical bridge index (e.g. site_idx_obs[i])
                  s_idx_var <- get_struct_index(s_name, response, hierarchical_info)

                  s_def <- jags_structure_definition(s_obj, variable_name = response, s_name = s_name, loop_bound = s_bound, zeros_name = s_zeros, category_index = "k", is_multi = is_struct_multi(s_name), i_index = s_idx_var)
                  if (!is.null(s_def)) {
                    model_lines <- safe_add_lines(model_lines, s_def$model_lines)
                    
                    if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
                       total_u <- paste0(total_u, " + ", s_def$term)
                    }

                    # Map structural parameters (DEDUPLICATED with Hierarchical Priority)
                    param_name <- paste0("sigma_", response, "_", s_name)
                    tau_name <- paste0("tau_u_", response, "_", s_name, "[k]")
                    
                    # Find if this signal pair already exists
                    existing_idx <- which(sapply(param_map, function(p) {
                      p$response == response && p$predictor == s_name && p$type == "structure" && grepl("\\[k\\]$", p$parameter)
                    }))
                    
                    # [UNIFICATION] Register the unified name if detected
                    is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(u, s_name), logical(1)))
                    
                    if (is_unified) {
                        unified_tau <- paste0("tau_u_", s_name, "_", response, "[k]")
                        unified_sigma <- paste0("sigma_", s_name, "_", response, "[k]")
                        param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                        param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
                    } else {
                        if (length(existing_idx) > 0) {
                          # [PRIORITY FIX] If existing is legacy and current is hierarchical, overwrite
                          existing_param <- param_map[[existing_idx[1]]]$parameter
                          if (grepl(paste0("_", s_name, "\\[k\\]$"), paste0(param_name, "[k]")) && !grepl(paste0("_", s_name, "\\[k\\]$"), existing_param)) {
                             param_map[[existing_idx[1]]]$parameter <- tau_name
                             param_map[[existing_idx[2]]]$parameter <- paste0(param_name, "[k]")
                          }
                        } else {
                          is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(u, s_name), logical(1)))
                          if (is_unified) {
                               unified_tau <- paste0("tau_u_", s_name, "_", response)
                               unified_sigma <- paste0("sigma_", s_name, "_", response)
                               param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                               param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0(unified_sigma, "[k]"), equation_index = NA, type = "structure")
                          } else {
                               param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = tau_name, equation_index = NA, type = "structure")
                               param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0(param_name, "[k]"), equation_index = NA, type = "structure")
                          }
                        }
                    }
                  }
                }
              }

              # 2. Random Group Structures
              for (r_name in random_structure_names) {
                if (!is_valid_random_level(response, r_name, hierarchical_info)) next
                
                s_suffix <- paste0("_", r_name)
                u_std <- paste0("u_std_", response, suffix, s_suffix)
                u <- paste0("u_", response, suffix, s_suffix)
                tau_u <- paste0("tau_u_", response, suffix, s_suffix)
                
                r_lvl <- get_random_level(response, r_name, hierarchical_info)
                # [SYNTAX GUARD] Final fallback to prevent empty [1:] indices
                n_groups <- if (!is.null(r_lvl) && nchar(r_lvl) > 0) paste0("N_", r_lvl) else if (nchar(r_name) > 0) paste0("N_", r_name) else "N"
                zeros_name <- if (!is.null(r_lvl) && nchar(r_lvl) > 0) paste0("zeros_", r_lvl) else if (nchar(r_name) > 0) paste0("zeros_", r_name) else "zeros"
                prec_name <- paste0("Prec_", r_name)

                model_lines <- safe_add_lines(model_lines, c(
                  paste0("  ", u_std, "[1:", n_groups, ", k] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
                  paste0("  for (g in 1:", n_groups, ") { ", u, "[g, k] <- ", u_std, "[g, k] / sqrt(", tau_u, "[k]) }")
                ))
                
                group_idx <- get_group_idx_string(response, r_name, hierarchical_info)
                total_u <- paste0(total_u, " + ", u, "[", group_idx, ", k]")
                
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0(tau_u, "[k]"), equation_index = NA, type = "structure")
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name, "[k]"), equation_index = NA, type = "structure")
              }

              # 3. Observation Loop for Multinomial (Inside category loop)
              model_lines <- c(
                model_lines,
                paste0("    for (i in 1:", loop_bound, ") {"),
                paste0("      ", epsilon, "[i, k] ~ dnorm(0, ", tau_res, "[k])"),
                paste0("      ", err, "[i, k] <- ", epsilon, "[i, k]", total_u),
                paste0("    }"),
                paste0("  }")
              )
            }
          } else if (dist == "ordinal") {
            # Ordinal error term: err[1:N]
            err <- paste0("err_", response, suffix)

            if (independent) {
              # Independent Ordinal: Standard GLM
              model_lines <- c(
                model_lines,
                paste0("  # Independent (Standard GLM) for ordinal: ", response),
                paste0("  for (i in 1:", loop_bound, ") {"),
                paste0("    ", err, "[i] <- 0"),
                paste0("  }")
              )
            } else {
              # Random Effects Formulation (Default/Optimised)
              epsilon <- paste0("epsilon_", response, suffix)
              tau_res <- paste0("tau_res_", response, suffix)
              
              # Initialize accumulators
              total_u <- ""

              # 1. Custom Structures (Hooks)
              if (length(structures) > 0) {
                processed_ord_signals <- character(0)
                for (s_idx in seq_along(structures)) {
                  s_name <- names(structures)[s_idx]
                  if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
                  
                  # [SINGULARITY GUARD]
                  if (s_name %in% processed_ord_signals) next
                  processed_ord_signals <- c(processed_ord_signals, s_name)
                  
                  s_obj <- structures[[s_idx]]

                  # [LINEAGE GUARD] Only process valid hierarchical branches
                  if (!is_valid_structure_mapping(get_struct_lvl(s_name, hierarchical_info), get_var_level(response, hierarchical_info), hierarchical_info, allow_identity = TRUE)) next

                  s_lvl <- get_struct_lvl(s_name, hierarchical_info)
                  s_bound <- if (is.null(s_lvl)) get_loop_bound(response, hierarchical_info) else paste0("N_", s_lvl)
                  s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)

                  # Hierarchical bridge index (e.g. site_idx_obs[i])
                  s_idx_var <- get_struct_index(s_name, response, hierarchical_info)

                  s_def <- jags_structure_definition(s_obj, variable_name = response, s_name = s_name, loop_bound = s_bound, zeros_name = s_zeros, is_multi = is_struct_multi(s_name), i_index = s_idx_var)
                  if (!is.null(s_def)) {
                    model_lines <- safe_add_lines(model_lines, s_def$model_lines)
                    
                    if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
                       total_u <- paste0(total_u, " + ", s_def$term)
                    }

                    # Map structural parameters (DEDUPLICATED)
                    param_name <- paste0("sigma_", response, "_", s_name)
                    tau_name <- paste0("tau_u_", response, "_", s_name)

                    is_duplicate <- any(sapply(param_map, function(p) {
                      p$response == response && p$predictor == s_name && p$type == "structure" && !grepl("\\[k\\]$", p$parameter)
                    }))

                    # [UNIFICATION] Register unified name in ordinal block
                    is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(u, s_name), logical(1)))
                    
                    if (is_unified) {
                        unified_tau <- paste0("tau_u_", s_name, "_", response)
                        unified_sigma <- paste0("sigma_", s_name, "_", response)
                        param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                        param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
                    } else {
                        if (length(existing_idx) > 0) {
                          # [PRIORITY FIX] If existing is legacy and current is hierarchical, overwrite
                          existing_param <- param_map[[existing_idx[1]]]$parameter
                          if (grepl(paste0("_", s_name, "$"), param_name) && !grepl(paste0("_", s_name, "$"), existing_param)) {
                             param_map[[existing_idx[1]]]$parameter <- tau_name
                             param_map[[existing_idx[2]]]$parameter <- param_name
                          }
                        } else {
                          param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = tau_name, equation_index = NA, type = "structure")
                          param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = param_name, equation_index = NA, type = "structure")
                        }
                    }
                  }
                }
              }

              # 2. Random Groups for Ordinal
              for (r_name in random_structure_names) {
                if (!is_valid_random_level(response, r_name, hierarchical_info)) next
                
                s_suffix <- paste0("_", r_name)
                u_std <- paste0("u_std_", response, suffix, s_suffix)
                u <- paste0("u_", response, suffix, s_suffix)
                tau_u <- paste0("tau_u_", response, suffix, s_suffix)
                
                # Determine hierarchical level if available
                r_lvl <- get_random_level(response, r_name, hierarchical_info)
                # [SYNTAX GUARD] Final fallback to prevent empty [1:] indices
                n_groups <- if (!is.null(r_lvl) && nchar(r_lvl) > 0) paste0("N_", r_lvl) else if (nchar(r_name) > 0) paste0("N_", r_name) else "N"
                zeros_name <- if (!is.null(r_lvl) && nchar(r_lvl) > 0) paste0("zeros_", r_lvl) else if (nchar(r_name) > 0) paste0("zeros_", r_name) else "zeros"
                prec_name <- paste0("Prec_", r_name)

                model_lines <- safe_add_lines(model_lines, c(
                  paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
                  paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- ", u_std, "[g] / sqrt(", tau_u, ") }")
                ))
                
                group_idx <- get_group_idx_string(response, r_name, hierarchical_info)
                total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
                
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = tau_u, equation_index = NA, type = "structure")
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name), equation_index = NA, type = "structure")
              }

              # 3. Observation Loop for Ordinal
              mag_terms <- vars_error_terms[[response]]
              total_mag <- if (length(mag_terms) > 0) paste0(" + ", paste(mag_terms, collapse = " + ")) else ""

              model_lines <- c(
                model_lines,
                paste0("  for (i in 1:", loop_bound, ") {"),
                paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_res, ")"),
                paste0("    ", err, "[i] <- ", epsilon, "[i]", total_u, total_mag),
                paste0("  }")
              )
            }
          } else if (dist == "poisson" || dist == "zip") {
        # Poisson error term: err[1:N]
        err <- paste0("err_", response, suffix)

        if (independent) {
          # Independent Poisson: Standard GLM (no overdispersion unless MAG present)
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) paste0(" + ", paste(mag_terms, collapse = " + ")) else ""

          model_lines <- c(
            model_lines,
            paste0("  # Independent (Standard GLM) for Poisson: ", response),
            paste0("  for (i in 1:", loop_bound, ") {"),
            paste0("    ", err, "[i] <- 0", total_mag),
            paste0("  }")
          )
        } else {
          # Optimized Random Effects
          epsilon <- paste0("epsilon_", response, suffix)
          tau_res <- paste0("tau_res_", response, suffix)
          
          # Initialize accumulators
          total_u <- ""
          local_obs_code <- c()

          # 1. Custom Structures (Hooks)
          if (length(structures) > 0) {
            for (s_idx in seq_along(structures)) {
              s_name <- names(structures)[s_idx]
              if (is.null(s_name) || s_name == "") s_name <- paste0("Struct", s_idx)
              s_obj <- structures[[s_idx]]

              # [LINEAGE GUARD] Only process valid hierarchical branches
              if (!is_valid_structure_mapping(get_struct_lvl(s_name, hierarchical_info), get_var_level(response, hierarchical_info), hierarchical_info, allow_identity = TRUE)) next

              s_lvl <- get_struct_lvl(s_name, hierarchical_info)
              s_bound <- if (is.null(s_lvl)) get_loop_bound(response, hierarchical_info) else paste0("N_", s_lvl)
              s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)

              # Hierarchical bridge index (e.g. site_idx_obs[i])
              s_idx_var <- get_struct_index(s_name, response, hierarchical_info)

              s_def <- jags_structure_definition(s_obj, variable_name = response, s_name = s_name, loop_bound = s_bound, zeros_name = s_zeros, is_multi = is_struct_multi(s_name), i_index = s_idx_var)
              if (!is.null(s_def)) {
                if (!is.null(s_def$setup_code)) model_lines <- safe_add_lines(model_lines, s_def$setup_code)
                
                # Global Prior Math
                model_lines <- safe_add_lines(model_lines, s_def$model_lines)
                
                # Observation Math Guard
                if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
                   total_u <- paste0(total_u, " + ", s_def$term)
                }

                # Map structural parameters
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
              }
            }
          }

          # 2. Random Group Structures
          for (r_name in random_structure_names) {
            # [LINEAGE GUARD] Only process valid hierarchical branches
            if (!is_valid_random_level(response, r_name, hierarchical_info)) next
            
            s_suffix <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)
            
            # Determine hierarchical level if available
            r_lvl <- get_random_level(response, r_name, hierarchical_info)
            n_groups <- if (!is.null(r_lvl)) paste0("N_", r_lvl) else paste0("N_", r_name)
            zeros_name <- if (!is.null(r_lvl)) paste0("zeros_", r_lvl) else paste0("zeros_", r_name)
            prec_name <- paste0("Prec_", r_name)

            model_lines <- c(
              model_lines,
              paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
              paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- ", u_std, "[g] / sqrt(", tau_u, ") }")
            )
            
            group_idx <- get_group_idx_string(response, r_name, hierarchical_info)
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
            
            param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = tau_u, equation_index = NA, type = "structure")
            param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name), equation_index = NA, type = "structure")
          }

          # 3. Final Observation Loop for Poisson
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) paste0(" + ", paste(mag_terms, collapse = " + ")) else ""

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for Poisson: ", response),
            paste0("  for (i in 1:", loop_bound, ") {"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_res, ")"),
            paste0("    ", err, "[i] <- ", epsilon, "[i]", total_u, total_mag),
            paste0("  }")
          )
        }
      } else if (dist == "occupancy") {
        # Occupancy Model (Single-Season)
        # Process:
        #   State: z[i] ~ dbern(psi[i]) where logit(psi) = mu_Response (reusing mu generation)
        #   Obs:   y[i,j] ~ dbern(z[i] * p[i,j]) where logit(p) = mu_p_Response

        # Determine if there is a detection model (p_Response)
        p_name <- paste0("p_", response)
        # Note: We assume mu_p_name is generated elsewhere if p_Response is in equations
        # JAGS is declarative, so order doesn't matter.

        # Initialize structure effect accumulator for psi
        total_u <- ""
        if (!is.null(structures)) {
          for (s_name in names(structures)) {
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) {
              if (!is.null(hierarchical_info)) {
                hierarchical_info$counts[[1]]
              } else {
                "N"
              }
            } else {
              paste0("N_", s_lvl)
            }

            # Call Generic
            s_def <- jags_structure_definition(
              structures[[s_name]],
              variable_name = response,
              s_name = s_name,
              loop_bound = s_bound,
              category_index = NULL
            )

            if (!is.null(s_def)) {
              model_lines <- safe_add_lines(model_lines, s_def$model_lines)
              total_u <- paste0(total_u, " + ", s_def$term)

              is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
              if (is_unified) {
                  unified_tau <- paste0("tau_u_", s_name, "_", response)
                  unified_sigma <- paste0("sigma_", s_name, "_", response)
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
              } else {
                  # Map structural parameters for post-processing/lambda
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
              }
            }
          }
        }

        # Start Observation Loop and Link Psi
        model_lines <- c(
          model_lines,
          paste0("  # Occupancy Model for ", response),
          # eq_loop_N is available here
          paste0("  for (i in 1:", eq_loop_N, ") {"),
          paste0(
            "    logit(psi_",
            response,
            "[i]) <- mu_",
            response,
            suffix,
            "[i]",
            total_u
          )
        )

        # Family Generic Dispatch
        fam_obj <- get_family_object(dist)

        # Locate predictors for this response (needed for generic)
        curr_pred <- NULL
        for (eq in eq_list) {
          if (eq$response == response) {
            curr_pred <- eq$predictors
            break
          }
        }

        def <- jags_family_definition(fam_obj, response, curr_pred)
        if (!is.null(def$model_code)) {
          model_lines <- safe_add_lines(model_lines, def$model_code)
        } else {
          stop(paste("Unknown distribution or missing module for:", dist))
        }

        model_lines <- c(model_lines, "  }")
      } else if (dist == "negbinomial" || dist == "zinb") {
        # Negative Binomial error term: err[1:N]
        # Single phylogenetic effect (like Poisson/ordinal)
        err <- paste0("err_", response, suffix)

        if (TRUE) {
          # Optimized Random Effects
          epsilon <- paste0("epsilon_", response, suffix)
          tau_res <- paste0("tau_res_", response, suffix)

          total_u <- ""

          # Custom Structures (Hooks)
          if (length(structures) > 0) {
            for (s_idx in seq_along(structures)) {
              s_name <- names(structures)[s_idx]
              if (is.null(s_name) || s_name == "") {
                s_name <- paste0("Struct", s_idx)
              }
              s_obj <- structures[[s_idx]]

              if (
                !is_valid_structure_mapping(
                  get_struct_lvl(s_name, hierarchical_info),
                  get_var_level(response, hierarchical_info),
                  hierarchical_info,
                  allow_identity = TRUE
                )
              ) {
                next
              }

              # Use Hook to get JAGS code bits
              s_def <- jags_structure_definition(
                s_obj,
                variable_name = response,
                s_name = s_name,
                loop_bound = get_loop_bound(response, hierarchical_info),
                is_multi = is_struct_multi(s_name)
              )

              if (!is.null(s_def)) {
                model_lines <- safe_add_lines(model_lines, s_def$model_lines)
                total_u <- paste0(total_u, " + ", s_def$term)

                # Map structural parameters for post-processing/lambda
                is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
                if (is_unified) {
                  unified_tau <- paste0("tau_u_", s_name, "_", response)
                  unified_sigma <- paste0("sigma_", s_name, "_", response)
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
                } else {
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("tau_u_", response, "_", s_name), equation_index = NA, type = "structure")
                  param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
                }
              }
            }
          }

          # Random Group Structures
          for (r_name in random_structure_names) {
            # [LINEAGE GUARD] Only process valid hierarchical branches
            if (!is_valid_random_level(response, r_name, hierarchical_info)) next
            
            s_suffix <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)
            
            r_lvl <- get_random_level(response, r_name, hierarchical_info)
            n_groups <- if (!is.null(r_lvl)) paste0("N_", r_lvl) else paste0("N_", r_name)
            zeros_name <- if (!is.null(r_lvl)) paste0("zeros_", r_lvl) else paste0("zeros_", r_name)
            prec_name <- paste0("Prec_", r_name)

            model_lines <- c(
              model_lines,
              paste0("  ", u_std, "[1:", n_groups, "] ~ dmnorm(", zeros_name, "[1:", n_groups, "], ", prec_name, "[1:", n_groups, ", 1:", n_groups, "])"),
              paste0("  for (g in 1:", n_groups, ") { ", u, "[g] <- ", u_std, "[g] / sqrt(", tau_u, ") }")
            )
            
            group_idx <- get_group_idx_string(response, r_name, hierarchical_info)
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
            
            param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = tau_u, equation_index = NA, type = "structure")
            param_map[[length(param_map) + 1]] <- list(response = response, predictor = r_name, parameter = paste0("sigma_", response, "_", r_name), equation_index = NA, type = "structure")
          }

          # For Negative Binomial / ZINB, we do NOT add an independent residual error (epsilon)
          # because the distribution already handles overdispersion via 'r' (size).
          # Adding epsilon creates identifiability issues (double overdispersion).

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for Negative Binomial: ", response),
            paste0("  for (i in 1:", eq_loop_N, ") {")
          )

          # Sum all terms: random effects + MAG terms
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) {
            paste0(" + ", paste(mag_terms, collapse = " + "))
          } else {
            ""
          }

          if (total_u == "" && total_mag == "") {
            model_lines <- safe_add_lines(model_lines, paste0("    ", err, "[i] <- 0"))
          } else {
            # total_u starts with " + ...", so we prepend "0"
            model_lines <- c(
              model_lines,
              paste0("    ", err, "[i] <- 0", total_u, total_mag)
            )
          }

          model_lines <- c(model_lines, "  }")
        }
      }
    }
  }

  # 2. Generate Likelihood for each variable (summing error terms)
  for (var in names(vars_error_terms)) {
    # Skip non-Gaussian variables (handled in GLMM error term blocks)
    dist <- if (!is.null(dist_list[[var]])) dist_list[[var]] else "gaussian"
    if (dist != "gaussian") {
      next
    }

    err_terms <- vars_error_terms[[var]]
    suffix <- if ((response_counter[[var]] %||% 0) > 1) "1" else ""
    loop_bound <- get_loop_bound(var, hierarchical_info)

    # Define Structure Error (Hooks)
    structure_term_str <- ""
    if (length(structures) > 0) {
      structure_terms <- character()
      for (s_idx in seq_along(structures)) {
        s_name <- names(structures)[s_idx]
        if (is.null(s_name) || s_name == "") {
          s_name <- paste0("Struct", s_idx)
        }
        s_obj <- structures[[s_idx]]

        # [PARTITIONING] Check if this variable qualifies for variance partitioning
        # (Exactly one structure + Gaussian + no other residual noise terms)
        local_struct_count <- 0
        for (sn in names(structures)) {
          if (is_valid_structure_mapping(get_struct_lvl(sn, hierarchical_info), get_var_level(var, hierarchical_info), hierarchical_info, allow_identity = TRUE)) {
            local_struct_count <- local_struct_count + 1
          }
        }
        
        # Partitioning triggers if Gaussian and exactly one structure is mapping
        use_partitioning_flag <- (dist == "gaussian" && local_struct_count == 1 && length(err_terms) == 0)

        # Use Hook to get JAGS code bits
        s_def <- jags_structure_definition(
          s_obj,
          variable_name = var,
          s_name = s_name,
          loop_bound = loop_bound,
          is_multi = is_struct_multi(s_name),
          use_partitioning = use_partitioning_flag
        )

        if (!is.null(s_def)) {
          model_lines <- safe_add_lines(model_lines, s_def$model_lines)
          structure_terms <- c(structure_terms, s_def$term)

          # Map structural parameters (DEDUPLICATED)
          param_name <- paste0("sigma_", var, "_", s_name)
          tau_name <- paste0("tau_u_", var, "_", s_name)

          is_duplicate <- any(sapply(param_map, function(p) {
            p$response == var && p$predictor == s_name && p$type == "structure"
          }))

          if (!is_duplicate) {
            param_map[[length(param_map) + 1]] <- list(
              response = var,
              predictor = s_name,
              parameter = tau_name,
              type = "structure"
            )
            param_map[[length(param_map) + 1]] <- list(
              response = var,
              predictor = s_name,
              parameter = param_name,
              type = "structure"
            )
          }
        }
      }
      if (length(structure_terms) > 0) {
        structure_term_str <- paste0(
          " + ",
          paste(structure_terms, collapse = " + ")
        )
      }
    }

    # --- UNIVERSAL VARIANCE PARTITIONING (Pagel's Lambda style) ---
    # Condition: Gaussian + Exactly one structure + No other noise terms (err_res)
    # This significantly improves identification by partitioning a single sigma_Total
    # instead of adding independent structural and residual variances.
    use_partitioning <- FALSE
    s_name_partition <- NULL
    if (dist == "gaussian" && length(structure_terms) == 1 && length(err_terms) == 0 && engine %in% c("jags", "nimble")) {
        # Identify the structure name (heuristic: find the first one in the loop)
        for (sn in names(structures)) {
            if (is_valid_structure_mapping(get_struct_lvl(sn, hierarchical_info), get_var_level(var, hierarchical_info), hierarchical_info, allow_identity = TRUE)) {
                s_name_partition <- sn
                use_partitioning <- TRUE
                break
            }
        }
    }

    if (use_partitioning && !is.null(s_name_partition)) {
        partition_param <- paste0("lambda_", var)
        sigma_total_param <- paste0("sigma_total_", var)
        tau_total_param <- paste0("tau_total_", var)
        
        # Override structural precision name to match partitioning logic
        tau_struct_partition <- paste0("tau_u_", var, "_", s_name_partition)
        tau_obs_partition <- paste0("tau_obs_", var)

        model_lines <- safe_add_lines(
            model_lines,
            c(
                paste0("  # Variance Partitioning (Pagel's Lambda style) for ", var),
                paste0("  ", sigma_total_param, " ~ dunif(0, 10)"),
                paste0("  ", partition_param, " ~ dunif(0, 1)"),
                paste0("  ", tau_total_param, " <- 1/(", sigma_total_param, " * ", sigma_total_param, ")"),
                paste0("  ", tau_struct_partition, " <- ", tau_total_param, " / max(0.001, ", partition_param, ")"),
                paste0("  ", tau_obs_partition, " <- ", tau_total_param, " / max(0.001, 1 - ", partition_param, ")")
            )
        )
        
        # Map partitioning parameters for post-processing
        param_map[[length(param_map) + 1]] <- list(response = var, predictor = s_name_partition, parameter = partition_param, type = "structure")
        param_map[[length(param_map) + 1]] <- list(response = var, predictor = s_name_partition, parameter = sigma_total_param, type = "structure")

    } else {
        # Define Observation Precision (Standard Additive Model)
        model_lines <- safe_add_lines(
            model_lines,
            c(
                get_precision_prior(paste0("tau_obs_", var), var),
                paste0("  sigma_obs_", var, " <- 1/sqrt(tau_obs_", var, ")")
            )
        )
    }

    # Build sum string
    sum_res_errs <- if (length(err_terms) > 0) {
        paste0(" + ", paste(err_terms, collapse = " + "))
    } else {
        ""
    }

    model_lines <- c(
      model_lines,
      paste0("  for (i in 1:", loop_bound, ") {"),
      paste0(
        "    ",
        var,
        "[i] ~ dnorm(mu_",
        var,
        suffix,
        "[i]",
        structure_term_str,
        sum_res_errs,
        ", tau_obs_",
        var,
        ")"
      )
    )
    if (engine == "jags") {
      model_lines <- c(
        model_lines,
        paste0(
          "    log_lik_",
          var,
          suffix,
          "[i] <- logdensity.norm(",
          var,
          "[i], mu_",
          var,
          suffix,
          "[i]",
          structure_term_str,
          " + ",
          sum_res_errs,
          ", tau_obs_",
          var,
          ")"
        )
      )
    }
    model_lines <- c(model_lines, paste0("  }"))
  }

  # Measurement error / Variability
  if (length(variability_list) > 0) {
    model_lines <- safe_add_lines(model_lines, "  # Measurement error / Variability")

    for (var in names(variability_list)) {
      type <- variability_list[[var]]

      # Skip measurement error block for occupancy variables (handled in likelihood)
      if (!is.null(dist_list[[var]]) && dist_list[[var]] == "occupancy") {
        next
      }

      # Ensure the variable is in the model
      if (!var %in% all_vars) {
        warning(paste(
          "Variable",
          var,
          "specified in 'variability' but not found in equations."
        ))
        next
      }

      if (type == "se") {
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:", eq_loop_N, ") {"),
          paste0(
            "    ",
            var,
            "_tau_obs[i] <- 1/(",
            var,
            "_se[i] * ",
            var,
            "_se[i])"
          ),
          paste0(
            "    ",
            var,
            "_mean[i] ~ dnorm(",
            var,
            "[i], ",
            var,
            "_tau_obs[i])"
          )
        )
        if (engine == "jags") {
          model_lines <- c(
            model_lines,
            paste0(
              "    log_lik_",
              var,
              "_mean[i] <- logdensity.norm(",
              var,
              "_mean[i], ",
              var,
              "[i], ",
              var,
              "_tau_obs[i])"
            )
          )
        }
        model_lines <- c(model_lines, paste0("  }"))
      } else if (type == "reps") {
        # For repeated measures, we need log_lik for EACH observation
        # Then sum them per individual for WAIC
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:", eq_loop_N, ") {")
        )
        # Replaced imperative accumulation with declarative summation
        model_lines <- c(
          model_lines,
          paste0(
            "    for (j in 1:",
            if (!is.null(hierarchical_info)) "N_reps_" else "N_reps_",
            var,
            "[i]) {"
          ),
          paste0(
            "      ",
            var,
            "_obs[i, j] ~ dnorm(",
            var,
            "[i], ",
            var,
            "_tau)"
          )
        )
        if (engine == "jags") {
          model_lines <- c(
            model_lines,
            paste0(
              "      lik_matrix_",
              var,
              "[i, j] <- logdensity.norm(",
              var,
              "_obs[i, j], ",
              var,
              "[i], ",
              var,
              "_tau)"
            )
          )
        }
        model_lines <- c(
          model_lines,
          paste0("    }")
        )
        if (engine == "jags") {
          model_lines <- c(
            model_lines,
            # Sum the log-likelihoods for this individual
            paste0(
              "    log_lik_",
              var,
              "_reps[i] <- sum(lik_matrix_",
              var,
              "[i, 1:N_reps_",
              var,
              "[i]])"
            )
          )
        }
        model_lines <- c(
          model_lines,
          paste0("  }")
        )
        model_lines <- safe_add_lines(
          model_lines,
          c(
            get_precision_prior(paste0(var, "_tau"), var),
            paste0("  ", var, "_sigma <- 1/sqrt(", var, "_tau)")
          )
        )
      } else {
        warning(paste("Unknown variability type:", type, "for variable", var))
      }
    }
  }

  model_lines <- safe_add_lines(model_lines, "  # Priors for structural parameters")

  # Priors for alpha, lambda, tau, sigma
  for (response in unique(names(response_counter))) {
    # Skip correlated vars (handled separately for Gaussian)
    dist <- dist_list[[response]] %||% "gaussian"
    if (response %in% correlated_vars && dist == "gaussian") {
      next
    }

    # Skip multinomial, ordinal (handled separately)
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "multinomial" || dist == "ordinal") {
      next
    }

    # Identify if this is an occupancy auxiliary parameter
    is_occupancy_aux <- FALSE
    if (grepl("^p_", response)) {
      target_resp <- sub("^p_", "", response)
      is_occupancy_aux <- any(sapply(names(dist_list), function(n) {
        (dist_list[[n]] == "occupancy") && (paste0("p_", n) == response)
      }))
    }

    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      # alpha prior
      alpha_name <- paste0("alpha_", response, suffix)
      # Determine default prior based on distribution (logit link needs tighter prior)
      logit_dists <- c(
        "occupancy",
        "binomial",
        "zip",
        "zinb",
        "bernoulli",
        "multinomial",
        "ordinal"
      )
      default_alpha <- "dnorm(0, 0.01)"
      if (dist %in% logit_dists || is_occupancy_aux) {
        default_alpha <- "dnorm(0, 1)"
      }
      model_lines <- c(
        model_lines,
        paste0("  ", get_prior(alpha_name, type = "alpha"))
      )

      # Only generate lambda/tau priors if NOT using GLMM (missing data),
      # UNLESS we are optimizing (which uses standard priors even for missing data)
      # Skip for occupancy models and their detection auxiliaries (no sigma/tau_e)
      # Condition to enter variance generation block
      # Basic rule: Skip if GLMM handles it (unless optimizing), or if distribution doesn't have residual variance.
      # HOWEVER, for standard GLMMs (occupancy, binomial) with structure, we STILL need to generate tau_u (random effect variance)
      # even if we don't need tau_e (residual).

      needs_residual_variance <- !dist %in%
        c(
          "occupancy",
          "binomial",
          "zip",
          "zinb",
          "poisson",
          "negbinomial",
          "bernoulli",
          "multinomial",
          "ordinal"
        )
      needs_random_variance <- !independent &&
        TRUE &&
        (length(structure_names) > 0 || length(random_structure_names) > 0) &&
        dist != "occupancy"

      if (
        (is.null(vars_with_na) || !response %in% vars_with_na || TRUE) &&
          (needs_residual_variance || needs_random_variance) &&
          !is_occupancy_aux
      ) {
        if (independent) {
          # Independent Priors (only tau_e)
          # We can also generate sigma for monitoring convenience
          if (
            !is.null(fix_residual_variance) &&
              (response %in%
                names(fix_residual_variance) ||
                length(fix_residual_variance) == 1)
          ) {
            val <- if (response %in% names(fix_residual_variance)) {
              fix_residual_variance[[response]]
            } else {
              fix_residual_variance[[1]]
            }
            prec <- 1 / val # Inverse variance
            model_lines <- c(
              model_lines,
              paste0(
                "  tau_res_",
                response,
                suffix,
                " <- ",
                prec,
                " # Fixed residual variance"
              )
            )
          } else {
            model_lines <- safe_add_lines(
              model_lines,
              paste0(
                "  ",
                get_precision_prior(
                  paste0("tau_res_", response, suffix),
                  response
                )
              )
            )
          }
        } else {
          # Component-wise: Estimate independent variance components
          # Note: We now use this even for single structure to match the restored model logic
          if (
            !is.null(fix_residual_variance) &&
              (response %in%
                names(fix_residual_variance) ||
                length(fix_residual_variance) == 1)
          ) {
            val <- if (response %in% names(fix_residual_variance)) {
              fix_residual_variance[[response]]
            } else {
              fix_residual_variance[[1]]
            }
            prec <- 1 / val # Inverse variance
            model_lines <- c(
              model_lines,
              paste0(
                "  tau_res_",
                response,
                suffix,
                " <- ",
                prec,
                " # Fixed residual variance"
              )
            )
          } else {
            # Only define tau_e if NOT Negative Binomial or ZINB
            # (since those use 'r' size parameter for dispersion)
            if (!dist %in% c("negbinomial", "zinb")) {
              model_lines <- safe_add_lines(
                model_lines,
                paste0(
                  "  ",
                  get_precision_prior(
                    paste0("tau_res_", response, suffix),
                    response
                  )
                )
              )
            }
          }

          # Generate Structure Priors
          processed_struct_signals <- character(0)
          for (s_name in structure_names) {
            if (s_name %in% processed_struct_signals) next
            processed_struct_signals <- c(processed_struct_signals, s_name)

            if (
              !is_valid_structure_mapping(
                get_struct_lvl(s_name, hierarchical_info),
                get_var_level(response, hierarchical_info),
                hierarchical_info,
                allow_identity = TRUE
              )
            ) {
              next
            }
            # Structure properties
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            # Always use structure name suffix to match model generation
            s_suffix <- paste0("_", s_name)

            is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
            if (is_unified) {
                tau_u <- paste0("tau_u_", s_name, "_", response)
                sigma_name <- paste0("sigma_", s_name, "_", response)
            } else {
                tau_u <- paste0("tau_u_", response, suffix, s_suffix)
                sigma_name <- paste0("sigma_", response, suffix, s_suffix)
            }

            model_lines <- safe_add_lines(
              model_lines,
              c(
                paste0("  ", get_precision_prior(tau_u, response)),
                paste0(
                  "  ", sigma_name,
                  " <- 1/sqrt(",
                  tau_u,
                  ")"
                )
              )
            )
          }

          # Generate lambda for compatibility if single structure
          # Only if tau_e exists (i.e. not NB/ZINB)
          if (
            length(structure_names) == 1 &&
              length(random_structure_names) == 0 &&
              !dist %in% c("negbinomial", "zinb")
          ) {
            s_name <- structure_names[1]
            s_suffix <- paste0("_", s_name)

            tau_u_name <- paste0("tau_u_", response, suffix, s_suffix)
            model_lines <- safe_add_lines(
              model_lines,
              paste0(
                "  lambda",
                response,
                suffix,
                " <- (1/",
                tau_u_name,
                ") / ((1/",
                tau_u_name,
                ") + (1/tau_res_",
                response,
                suffix,
                "))"
              )
            )
          }

          # Random Effects Priors
          processed_rand_signals <- character(0)
          for (r_name in random_structure_names) {
            if (r_name %in% processed_rand_signals) next
            processed_rand_signals <- c(processed_rand_signals, r_name)

            # [LINEAGE GUARD] Only generate priors for valid hierarchical branches
            if (!is_valid_random_level(response, r_name, hierarchical_info)) next
            
            s_suffix <- paste0("_", r_name)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            model_lines <- safe_add_lines(
              model_lines,
              c(
                paste0("  ", get_precision_prior(tau_u, response)),
                paste0(
                  "  sigma_",
                  response,
                  suffix,
                  "_",
                  r_name,
                  " <- 1/sqrt(",
                  tau_u,
                  ")"
                )
              )
            )
          }

        }
      }
    }
  }

  # Priors for multinomial parameters (arrays)
  for (response in unique(names(response_counter))) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "multinomial") {
      K_var <- paste0("K_", response)
      if (independent) {
        if (
          !is.null(fix_residual_variance) &&
            (response %in%
              names(fix_residual_variance) ||
              length(fix_residual_variance) == 1)
        ) {
          val <- if (response %in% names(fix_residual_variance)) {
            fix_residual_variance[[response]]
          } else {
            fix_residual_variance[[1]]
          }
          prec <- 1 / val
          tau_line <- paste0(
            "    tau_res_",
            response,
            "[k] <- ",
            prec,
            " # Fixed"
          )
        } else {
          tau_line <- paste0("    ", get_precision_prior(paste0("tau_res_", response, "[k]"), response))
        }

        model_lines <- c(
          model_lines,
          paste0("  # Independent Priors for ", response, " (Multinomial)"),
          paste0("  alpha_", response, "[1] <- 0"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 0.01)"),
          tau_line,
          "  }"
        )
      } else {
        if (
          !is.null(fix_residual_variance) &&
            (response %in%
              names(fix_residual_variance) ||
              length(fix_residual_variance) == 1)
        ) {
          val <- if (response %in% names(fix_residual_variance)) {
            fix_residual_variance[[response]]
          } else {
            fix_residual_variance[[1]]
          }
          prec <- 1 / val
          tau_line <- paste0(
            "    tau_res_",
            response,
            "[k] <- ",
            prec,
            " # Fixed"
          )
        } else {
          tau_line <- paste0("    ", get_precision_prior(paste0("tau_res_", response, "[k]"), response))
        }

        model_lines <- c(
          model_lines,
          # Priors for ", response, " (Multinomial)"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 0.01)"),

          # Residual error
          tau_line,

          # Random effects priors (loop over structures)
          {
            prior_lines <- c()
            processed_multi_struct <- character(0)

            # Phylogenetic / N-dim Structures
            for (s_name in structure_names) {
              if (s_name %in% processed_multi_struct) next
              processed_multi_struct <- c(processed_multi_struct, s_name)

              if (
                !is_valid_structure_mapping(
                  get_struct_lvl(s_name, hierarchical_info),
                  get_var_level(response, hierarchical_info),
                  hierarchical_info,
                  allow_identity = TRUE
                )
              ) {
                next
              }
              # Structure properties
              s_lvl <- get_struct_lvl(s_name, hierarchical_info)
              s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
              s_zeros <- if (is.null(s_lvl)) {
                "zeros"
              } else {
                paste0("zeros_", s_lvl)
              }
              s_suffix <- paste0("_", s_name)
              tau_u <- paste0("tau_u_", response, s_suffix)

              prior_lines <- c(
                prior_lines,
                paste0(
                  "    ",
                  get_precision_prior(paste0(tau_u, "[k]"), response)
                ),
                paste0(
                  "    sigma_",
                  response,
                  s_suffix,
                  "[k] <- 1/sqrt(",
                  tau_u,
                  "[k])"
                )
              )
            }

            # Random Group Structures
            processed_multi_rand <- character(0)
            for (r_name in random_structure_names) {
              if (r_name %in% processed_multi_rand) next
              processed_multi_rand <- c(processed_multi_rand, r_name)

              s_suffix <- paste0("_", r_name)
              tau_u <- paste0("tau_u_", response, s_suffix)

              prior_lines <- c(
                prior_lines,
                paste0(
                  "    ",
                  get_precision_prior(paste0(tau_u, "[k]"), response)
                ),
                paste0(
                  "    sigma_",
                  response,
                  s_suffix,
                  "[k] <- 1/sqrt(",
                  tau_u,
                  "[k])"
                )
              )
            }
            prior_lines
          },

          # Derived lambda (only if single structure, for backward compatibility or convenience)
          if (
            length(structure_names) == 1 && length(random_structure_names) == 0
          ) {
            s_name <- structure_names[1]
            s_suffix <- paste0("_", s_name)
            tau_u_name <- paste0("tau_u_", response, s_suffix)
            paste0(
              "    lambda_",
              response,
              "[k] <- (1/",
              tau_u_name,
              "[k]) / ((1/",
              tau_u_name,
              "[k]) + (1/tau_res_",
              response,
              "[k]))"
            )
          }
        )
        model_lines <- c(
          model_lines,
          paste0("  # Priors for ", response, " (Multinomial)"),
          paste0("  alpha_", response, "[1] <- 0"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 0.01)"),
          paste0("    lambda_", response, "[k] ~ dunif(0, 1)"),
          paste0("    ", get_precision_prior(paste0("tau_", response, "[k]"), response)),
          "  }"
        )
      }

      # Betas (arrays)
      for (eq in eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            beta_name <- paste0("beta_", response, "_", pred)
            model_lines <- c(
              model_lines,
              paste0("  ", beta_name, "[1] <- 0"),
              paste0("  for (k in 2:", K_var, ") {"),
              paste0(
                "    ",
                get_prior(paste0(beta_name, "[k]"), type = "beta")
              ),
              "  }"
            )
          }
        }
      }
    }
  }

  # Priors for ordinal parameters (cutpoints + variance components)
  for (response in unique(names(response_counter))) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "ordinal") {
      K_var <- paste0("K_", response)

      # Loop over response instances (if there are repeats)
      for (k in 1:response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)

        model_lines <- c(
          model_lines,
          paste0("  # Priors for ", response, suffix, " (Ordinal)"),
          # Ordered cutpoints using delta transformation
          paste0("  cutpoint_raw_", response, suffix, "[1] ~ dnorm(0, 0.1)"),
          paste0(
            "  cutpoint_",
            response,
            suffix,
            "[1] <- cutpoint_raw_",
            response,
            suffix,
            "[1]"
          ),
          paste0("  for (k in 2:(", K_var, "-1)) {"),
          paste0("    cutpoint_raw_", response, suffix, "[k] ~ dnorm(0, 0.1)"),
          paste0(
            "    cutpoint_",
            response,
            suffix,
            "[k] <- cutpoint_",
            response,
            suffix,
            "[k-1] + exp(cutpoint_raw_",
            response,
            suffix,
            "[k])"
          ),
          "  }",
          # Variance components
          if (independent) {
            # Independent: Only tau_e
            c(
              paste0(
                "  ",
                get_precision_prior(
                  paste0("tau_res_", response, suffix),
                  response
                )
              ),
              paste0("  # No tau_u or lambda for independent model")
            )
          } else {
            # Structured: tau_u and tau_e and lambda
            c(
              paste0(
                "  ",
                get_precision_prior(
                  paste0("tau_u_", response, suffix),
                  response
                )
              ),
              paste0(
                "  ",
                get_precision_prior(
                  paste0("tau_res_", response, suffix),
                  response
                )
              ),
              # Derived lambda
              paste0(
                "  lambda_",
                response,
                suffix,
                " <- (1/tau_u_",
                response,
                suffix,
                ") / ((1/tau_u_",
                response,
                suffix,
                ") + (1/tau_res_",
                response,
                suffix,
                "))"
              )
            )
          }
        )
      }

      # Betas for ordinal predictors
      for (eq in eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            beta_name <- paste0("beta_", response, "_", pred)
            model_lines <- c(
              model_lines,
              paste0("  ", get_prior(beta_name, type = "beta"))
            )
          }
        }
      }
    }
  }

  # Priors for Negative Binomial parameters (size only - others handled in main loop)
  for (response in unique(names(response_counter))) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist == "negbinomial" || dist == "zinb") {
      # Loop over response instances
      for (k in 1:response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)

        model_lines <- c(
          model_lines,
          paste0(
            "  # Priors for ",
            response,
            suffix,
            " (Negative Binomial Size)"
          ),
          # Size parameter (controls overdispersion)
          paste0(
            "  r_",
            response,
            suffix,
            " ~ dgamma(0.01, 0.01)  # Vague prior for size"
          )
        )
      }
    }
  }

  # Priors for correlated vars alphas (intercepts - handled here for Gaussian only)
  for (var in correlated_vars) {
    # Skip non-Gaussian (handled in main priors loop)
    dist <- dist_list[[var]] %||% "gaussian"
    if (dist != "gaussian") {
      next
    }

    model_lines <- c(
      model_lines,
      paste0("  alpha_", var, " ~ dnorm(0, 0.01)")
    )
  }

  # Priors for Zero-Inflation parameters (psi)
  for (response in unique(names(response_counter))) {
    dist <- dist_list[[response]] %||% "gaussian"
    if (dist %in% c("zip", "zinb")) {
      for (k in 1:response_counter[[response]]) {
        suffix <- if (k == 1) "" else as.character(k)
        model_lines <- c(
          model_lines,
          paste0(
            "  psi_",
            response,
            suffix,
            " ~ dunif(0, 1) # Zero-inflation probability"
          )
        )
      }
    }
  }

  # Priors for regression coefficients
  unique_betas <- unique(unlist(beta_counter))

  # Exclude betas that are defined in distribution-specific sections
  # (multinomial, ordinal, poisson, negbinomial all define their own betas)
  excluded_betas <- c()

  for (response in unique(names(response_counter))) {
    dist <- dist_list[[response]] %||% "gaussian"

    # Skip likelihood generation for detection probability equations (p_Response)
    # These are auxiliary equations for Occupancy models
    if (grepl("^p_", response)) {
      target_resp <- sub("^p_", "", response)
      # Check if this target is an occupancy model response
      is_occupancy_aux <- any(sapply(names(dist_list), function(n) {
        (dist_list[[n]] == "occupancy") && (paste0("p_", n) == response)
      }))

      if (is_occupancy_aux) {
        # This is an auxiliary equation for detection probability
        # We generate the linear predictor (mu_p_Response) but NOT a likelihood
        # The likelihood is handled in the main occupancy block
        next
      }
    }

    if (dist %in% c("multinomial", "ordinal")) {
      for (eq in eq_list) {
        if (eq$response == response) {
          for (pred in eq$predictors) {
            excluded_betas <- c(
              excluded_betas,
              paste0("beta_", response, "_", pred)
            )
          }
        }
      }
    }
  }

  unique_betas <- setdiff(unique_betas, excluded_betas)

  # Track which latent variables have been "pinned" with a sign constraint on their first indicator
  pinned_latents <- character()

  for (beta in unique_betas) {
    # Isolate response part of beta name: beta_Response_Predictor
    # Assuming names are strictly beta_Resp_Pred, but Resp could contain underscores if user did so.
    # A safer way is using param_map, but we are in a simple loop here.
    # Let's try to match against names(dist_list).

    default_beta <- "dnorm(0, 0.01)"

    # Simple heuristic to find response name in beta string
    found_resp <- NULL
    for (nm in names(dist_list)) {
      # Check if beta starts with "beta_nm_"
      if (startsWith(beta, paste0("beta_", nm, "_"))) {
        # Valid match, but beware of prefixes (e.g. beta_vars vs beta_var)
        # Sort names by length descending to match longest possible response name first?
        # For now, let's assume exact match plus underscore.
        found_resp <- nm
        break
      }
      # Also check "beta_p_nm_" for occupancy auxiliary
      if (startsWith(beta, paste0("beta_p_", nm, "_"))) {
        # It's an occupancy detection parameter
        # These are logit scale -> need tight priors
        default_beta <- "dnorm(0, 1)"
        break
      }
    }

    if (!is.null(found_resp)) {
      d <- dist_list[[found_resp]] %||% "gaussian"
      if (
        d %in%
          c(
            "occupancy",
            "binomial",
            "zip",
            "zinb",
            "bernoulli",
            "multinomial",
            "ordinal"
          )
      ) {
        default_beta <- "dnorm(0, 1)"
      }
    }

    # APPLY ANCHORING FOR LATENT IDENTIFIABILITY
    # We use Unit Loading Identification (ULI) by default (fix first loading to 1.0).
    # This provides the most stable identification for non-linear GLVMs (Poisson/Binomial).
    if (!is.null(latent) && length(latent) > 0) {
      for (lat in latent) {
        # Match beta_ANYTHING_LatentVariable
        if (!lat %in% pinned_latents && grepl(paste0("_", lat, "$"), beta)) {
          if (fix_latent == "loading") {
            default_beta <- "1.0" # Unit Loading Identification (ULI)
          } else if (fix_latent == "sign") {
            default_beta <- paste0(default_beta, " T(0,)") # Sign Identification
          }
          pinned_latents <- c(pinned_latents, lat)
          break
        }
      }
    }

    model_lines <- c(
      model_lines,
      paste0("  ", get_prior(beta, type = "beta", default = default_beta))
    )
  }

  # Emit any extra prior declarations from extension packages (e.g., community
  # hyperparameters from because.occupancy's build_community_hyperpriors()).
  # Convention: entries in `priors` whose name starts with "__community__"
  # are emitted verbatim as JAGS lines (the value is a raw JAGS statement).
  if (!is.null(priors)) {
    community_extras <- priors[startsWith(names(priors), "__community__")]
    if (length(community_extras) > 0) {
      model_lines <- c(
        model_lines,
        "  # Community hyperparameter declarations (from extension package)"
      )
      for (line in unlist(community_extras)) {
        model_lines <- safe_add_lines(model_lines, paste0("  ", line))
      }
    }
  }

  # Object selection (if is_multi_structure)
  if (is_multi_structure) {
    model_lines <- c(
      model_lines,
      "  # Phylogenetic uncertainty weighting",
      "  for (k in 1:Ntree) {",
      "    p_tree[k] <- 1/Ntree",
      "  }",
      "  K ~ dcat(p_tree[1:Ntree])"
    )
  }

  # Covariance structures for responses

  for (response in names(response_counter)) {
    # Skip correlated vars (handled separately for Gaussian)
    dist <- dist_list[[response]] %||% "gaussian"
    if (response %in% correlated_vars && dist == "gaussian") {
      next
    }

    for (k in 1:response_counter[[response]]) {
      suffix <- if (k == 1) "" else as.character(k)
      dist <- dist_list[[response]] %||% "gaussian"

      # Skip TAU matrix for variables with missing data (using element-wise) or binomial error terms
      # Handle legacy GLMM blocking
      use_glmm <- (!is.null(vars_with_na) &&
        response %in% vars_with_na &&
        FALSE)

      if (dist == "gaussian" && !use_glmm) {
        if (FALSE) {
          if (is_multi_structure) {
            model_lines <- c(
              model_lines,
              paste0(
                "  Mlam",
                response,
                suffix,
                " <- lambda",
                response,
                suffix,
                "*multiVCV[,,K] + (1-lambda",
                response,
                suffix,
                ")*ID"
              ),
              paste0(
                "  TAU",
                tolower(response),
                suffix,
                " <- tau",
                response,
                suffix,
                "*inverse(Mlam",
                response,
                suffix,
                ")"
              )
            )
          } else {
            model_lines <- c(
              model_lines,
              paste0(
                "  Mlam",
                response,
                suffix,
                " <- lambda",
                response,
                suffix,
                "*VCV + (1-lambda",
                response,
                suffix,
                ")*ID"
              ),
              paste0(
                "  TAU",
                tolower(response),
                suffix,
                " <- tau",
                response,
                suffix,
                "*inverse(Mlam",
                response,
                suffix,
                ")"
              )
            )
          }
        }
      } else if (dist == "gaussian" && use_glmm) {
        # GLMM covariance for latent error term with generic structure support
        # err ~ dmnorm(0, tau_<structure> * inv(VCV))

        s_name <- if (length(structure_names) > 0) {
          structure_names[1]
        } else {
          "struct"
        }

        err <- paste0("err_", response, suffix)
        mu_err <- paste0("mu_err_", response, suffix)
        # [UNIFICATION] Silence legacy structural registration if dispatch already handled it
        tau_struct <- paste0("tau_u_", response, suffix, "_", s_name)
        tau_res <- paste0("tau_res_", response, suffix)

        # Only add legacy prior if not a unified structure
        is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(u, s_name), logical(1)))
        
        # [PARTITIONING] Check if this qualifies for variance partitioning
        # Exactly one structure + Gaussian + no other residual noise terms
        # (Heuristic: local_struct_count == 1)
        local_struct_count <- 0
        if (!is.null(structures)) {
          for (sn in names(structures)) {
            if (is_valid_structure_mapping(get_struct_lvl(sn, hierarchical_info), get_var_level(response, hierarchical_info), hierarchical_info, allow_identity = TRUE)) {
              local_struct_count <- local_struct_count + 1
            }
          }
        }
        use_partitioning <- (local_struct_count == 1 && engine %in% c("jags", "nimble"))

        if (use_partitioning) {
            partition_param <- paste0("lambda_", response)
            sigma_total_param <- paste0("sigma_total_", response)
            tau_total_param <- paste0("tau_total_", response)
            
            # Use the already defined prefix names for the partition components
            # (Ensures compatibility with s_def term generation)
            model_lines <- safe_add_lines(
                model_lines,
                c(
                    paste0("  # GLMM Variance Partitioning for ", response),
                    paste0("  ", sigma_total_param, " ~ dunif(0, 10)"),
                    paste0("  ", partition_param, " ~ dunif(0, 1)"),
                    paste0("  ", tau_total_param, " <- 1/(", sigma_total_param, " * ", sigma_total_param, ")"),
                    # Precision for structured effect (tau_struct)
                    paste0("  ", tau_struct, " <- ", tau_total_param, " / max(0.001, ", partition_param, ")"),
                    # Precision for residual error (tau_res)
                    paste0("  ", tau_res, " <- ", tau_total_param, " / max(0.001, 1 - ", partition_param, ")")
                )
            )
            # Map partitioning parameters
            param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = partition_param, type = "structure")
            param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = sigma_total_param, type = "structure")
            
        } else {
            # Standard Additive Priors
            model_lines <- safe_add_lines(
              model_lines,
              paste0("  ", get_precision_prior(tau_res, response))
            )

            if (is_unified) {
                # [UNIFICATION] Map the unified name for summary reporting
                unified_tau <- paste0("tau_u_", s_name, "_", response)
                unified_sigma <- paste0("sigma_", s_name, "_", response)
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_tau, equation_index = NA, type = "structure")
                param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = unified_sigma, equation_index = NA, type = "structure")
            } else {
                model_lines <- c(
                    model_lines,
                    paste0("  ", tau_struct, " ~ dgamma(10, 10)"),
                    # Calculate lambda for reporting
                    paste0(
                        "  lambda_",
                        response,
                        suffix,
                        "_",
                        s_name,
                        " <- (1/",
                        tau_struct,
                        ") / ((1/",
                        tau_struct,
                        ") + (1/",
                        tau_res,
                        "))"
                    )
                )
                if (!is_unified) {
                    param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = tau_struct, equation_index = NA, type = "structure")
                    param_map[[length(param_map) + 1]] <- list(response = response, predictor = s_name, parameter = paste0("sigma_", response, "_", s_name), equation_index = NA, type = "structure")
                }
            }
        }

        if (is_struct_multi(s_name)) {
          model_lines <- c(
            model_lines,
            paste0(
              "  TAU_",
              response,
              suffix,
              "_",
              s_name,
              " <- ",
              tau_struct,
              " * inverse(multiVCV[,,K])"
            ),
            paste0(
              "  ",
              err,
              "[1:loop_bound] ~ dmnorm(",
              mu_err,
              "[], TAU_",
              response,
              suffix,
              "_",
              s_name,
              ")"
            )
          )
        } else {
          model_lines <- c(
            model_lines,
            paste0(
              "  TAU_",
              response,
              suffix,
              "_",
              s_name,
              " <- ",
              tau_struct,
              " * inverse(VCV)"
            ),
            paste0(
              "  ",
              err,
              "[1:loop_bound] ~ dmnorm(",
              mu_err,
              "[], TAU_",
              response,
              suffix,
              "_",
              s_name,
              ")"
            )
          )
        }
      } else if (dist == "binomial") {
        # GLMM covariance for binomial error term
        if (FALSE) {
          # Marginal Formulation
          if (is_multi_structure) {
            model_lines <- c(
              model_lines,
              paste0(
                "  Mlam",
                response,
                suffix,
                " <- lambda",
                response,
                suffix,
                "*multiVCV[,,K] + (1-lambda",
                response,
                suffix,
                ")*ID"
              ),
              paste0(
                "  TAU",
                tolower(response),
                suffix,
                " <- tau",
                response,
                suffix,
                "*inverse(Mlam",
                response,
                suffix,
                ")"
              )
            )
          } else {
            model_lines <- c(
              model_lines,
              paste0(
                "  Mlam",
                response,
                suffix,
                " <- lambda",
                response,
                suffix,
                "*VCV + (1-lambda",
                response,
                suffix,
                ")*ID"
              ),
              paste0(
                "  TAU",
                tolower(response),
                suffix,
                " <- tau",
                response,
                suffix,
                "*inverse(Mlam",
                response,
                suffix,
                ")"
              )
            )
          }
        }
        # tau_u priors for optimized model are already generated in the main priors block
        # No need to duplicate here.
      } else if (dist == "occupancy") {
        # Occupancy priors are handled by the generic loop above (lines 2243+ for structures)
        # No specific residual variance (tau_e) needed.
      } else if (dist == "multinomial") {
        # Multinomial covariance
        # We need TAU[,,k] for each k
        K_var <- paste0("K_", response)

        if (FALSE) {
          # k=1 is reference category (fixed to identity)
          # k>=2 have estimated phylogenetic signal
          if (is_multi_structure) {
            model_lines <- c(
              model_lines,
              paste0("  # Covariance matrices for multinomial"),
              "  # Reference category k=1",
              paste0("  for (i in 1:", loop_bound, ") {"),
              paste0("    for (j in 1:", loop_bound, ") {"),
              paste0("      Mlam_", response, "[i,j,1] <- ID[i,j]"),
              "    }",
              "  }",
              paste0(
                paste0(
                  "  TAU_",
                  tolower(response),
                  "_",
                  suffix,
                  "[1:",
                  loop_bound,
                  ",1:",
                  loop_bound,
                  ",1] <- ID[1:",
                  loop_bound,
                  ",1:",
                  loop_bound,
                  "]"
                )
              ),
              "  # Estimated categories k>=2",
              paste0("  for (k in 2:", K_var, ") {"),
              paste0("    for (i in 1:", loop_bound, ") {"),
              paste0("      for (j in 1:", loop_bound, ") {"),
              paste0(
                "        Mlam_",
                response,
                "[i,j,k] <- lambda_",
                response,
                "[k]*multiVCV[i,j,K] + (1-lambda_",
                response,
                "[k])*ID[i,j]"
              ),
              "      }",
              "    }",
              paste0(
                "    TAU_",
                tolower(response),
                "_",
                suffix,
                "[1:",
                loop_bound,
                ",1:",
                loop_bound,
                ",k] <- inverse(Mlam_",
                response,
                "[,,k])"
              ),
              "  }"
            )
          }
        }
      }
    }
  }

  # Covariance for correlated vars (phylogenetic part)
  # Optimised models use eigendecomposition
  if (!is.null(induced_correlations) && FALSE) {
    s_name <- if (length(structure_names) > 0) structure_names[1] else "struct"
    for (var in correlated_vars) {
      if (is_struct_multi(s_name)) {
        model_lines <- c(
          model_lines,
          paste0(
            "  TAU_",
            s_name,
            "_",
            var,
            " <- tau_",
            s_name,
            "_",
            var,
            " * inverse(multiVCV[,,K])"
          )
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0(
            "  TAU_",
            s_name,
            "_",
            var,
            " <- tau_",
            s_name,
            "_",
            var,
            " * inverse(VCV)"
          )
        )
      }
    }
  }

  # Imputation priors for predictors (those not modeled as responses)
  if (
    length(setdiff(all_vars, names(response_counter))) > 0 &&
      (length(latent) > 0 ||
        (!is.null(vars_with_na) && length(vars_with_na) > 0))
  ) {
    model_lines <- safe_add_lines(model_lines, "  # Predictor priors for imputation")
  }
  non_response_vars <- setdiff(all_vars, names(response_counter))
  for (var in unique(non_response_vars)) {
    # Check if this is a latent variable
    is_latent <- !is.null(latent) && var %in% latent

    # Skip fully observed predictors (not latent and no missing data)
    # Only generate imputation for: latent variables OR variables with missing data
    if (!is_latent && (is.null(vars_with_na) || !var %in% vars_with_na)) {
      next # Skip this predictor - it's fully observed data
    }

    model_lines <- safe_add_lines(
      model_lines,
      c(
        paste0("  for (i in 1:", get_loop_bound(var, hierarchical_info), ") {"),
        paste0("    mu", var, "[i] <- 0"),
        paste0("  }")
      )
    )

    if (independent) {
      # Independent imputation (i.i.d normal)
      if (is_latent && standardize_latent) {
        # Use N(0,1) prior for standardized latent variable
        model_lines <- safe_add_lines(
          model_lines,
          c(
            paste0(
              "  for (i in 1:",
              get_loop_bound(var, hierarchical_info),
              ") {"
            ),
            paste0(
              "    ",
              var,
              "[i] ~ dnorm(mu", var, "[i], 1)  # Standardized latent variable with parents"
            ),
            paste0("  }")
          )
        )
      } else {
        model_lines <- safe_add_lines(
          model_lines,
          c(
            paste0(
              "  for (i in 1:",
              get_loop_bound(var, hierarchical_info),
              ") {"
            ),
            paste0("    ", var, "[i] ~ dnorm(mu", var, "[i], tau_res_", var, ")"),
            paste0("  }")
          )
        )
      }
    } else {
      # Optimized Random Effects Formulation for Predictors (Additive)

      # If latent variable with standardize_latent = TRUE, use simple N(0,1) prior
      # ignoring structure (assumes latent is standardized white noise)
      if (is_latent && standardize_latent) {
        model_lines <- c(
          model_lines,
          paste0(
            "  for (i in 1:",
            get_loop_bound(var, hierarchical_info),
            ") {"
          ),
          paste0(
            "    ",
            var,
            "[i] ~ dnorm(mu", var, "[i], 1)  # Standardized latent variable with parents"
          ),
          paste0("  }")
        )
      } else {
        # Standard random effects formulation
        additive_terms <- ""

        if (!is.null(structures)) {
          for (s_name in names(structures)) {
            # [LINEAGE GUARD] Skip structures that do not apply to this variable's level
            if (!is_valid_structure_mapping(get_struct_lvl(s_name, hierarchical_info), get_var_level(var, hierarchical_info), hierarchical_info, allow_identity = TRUE)) {
              next
            }
            
            s_suffix <- paste0("_", s_name)
            u_var <- paste0("u_", var, s_suffix)
            tau_u <- paste0("tau_u_", var, s_suffix)

            # Hierarchical info
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) get_loop_bound(var, hierarchical_info) else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)

            # Hierarchical bridge index (e.g. site_idx_obs[i])
            s_idx_var <- get_struct_index(s_name, var, hierarchical_info)

            # Determine if we should use partitioning (Lambda logic)
            r_lvl <- get_var_level(var, hierarchical_info)
            local_struct_count <- 0
            for (sn in names(structures)) {
              if (is_valid_structure_mapping(get_struct_lvl(sn, hierarchical_info), r_lvl, hierarchical_info, allow_identity = TRUE)) {
                local_struct_count <- local_struct_count + 1
              }
            }
            use_partitioning <- (local_struct_count == 1) && 
                               (get_family_object(dist_list[[var]] %||% "gaussian")$family == "gaussian") &&
                               (s_lvl == r_lvl)

            # Call Generic to define the error vector u ~ dmnorm(...)
            s_def <- jags_structure_definition(
              structures[[s_name]],
              variable_name = var,
              s_name = s_name,
              loop_bound = s_bound,
              zeros_name = s_zeros,
              is_multi = is_struct_multi(s_name),
              i_index = s_idx_var,
              use_partitioning = use_partitioning
            )

            if (!is.null(s_def$model_lines)) {
              model_lines <- safe_add_lines(model_lines, s_def$model_lines)
            }

            if (length(s_def$term) > 0 && nchar(s_def$term) > 0) {
               additive_terms <- paste0(additive_terms, " + ", s_def$term)
            }
          }
        }

        tau_e <- paste0("tau_res_", var)

        model_lines <- c(
          model_lines,
          paste0(
            "  for (i in 1:",
            get_loop_bound(var, hierarchical_info),
            ") {"
          ),
          paste0(
            "    ",
            var,
            "[i] ~ dnorm(mu",
            var,
            "[i]",
            additive_terms,
            ", ",
            tau_e,
            ")"
          ),
          paste0("  }")
        )
      }
    }

    # For latent variables, fix tau = 1 (standardize)
    # For observed predictors, estimate tau
    # Skip if standardize_latent = TRUE (variance already set in N(0,1) prior)
    if (is_latent && !standardize_latent) {
      model_lines <- c(
        model_lines,
        paste0("  # Latent variable: standardized (var = 1)"),
        paste0("  lambda", var, " ~ dunif(0, 1)"),
        paste0("  tau_u_", var, " <- 1/lambda", var),
        paste0("  tau_res_", var, " <- 1/(1-lambda", var, ")"),
        paste0("  sigma", var, " <- 1") # Fixed to 1
      )
    } else if (is_latent && standardize_latent) {
      # No variance parameters needed - already specified in N(0,1) prior
      model_lines <- c(
        model_lines,
        paste0(
          "  # Latent variable: fully standardized with N(0,1) prior (no variance parameters)"
        )
      )
    } else {
      if (independent) {
        # --- UNIVERSAL VARIANCE PARTITIONING (Independent Variable) ---
        # Check if we have exactly one structure active for this exogenous var
        local_struct_count <- 0
        s_name_partition <- NULL
        if (!is.null(structures)) {
            for (sn in names(structures)) {
                if (is_valid_structure_mapping(get_struct_lvl(sn, hierarchical_info), get_var_level(var, hierarchical_info), hierarchical_info, allow_identity = TRUE)) {
                    local_struct_count <- local_struct_count + 1
                    s_name_partition <- sn
                }
            }
        }

        if (local_struct_count == 1 && !is.null(s_name_partition)) {
            partition_param <- paste0("lambda_", var)
            sigma_total_param <- paste0("sigma_total_", var)
            tau_total_param <- paste0("tau_total_", var)
            
            tau_struct_partition <- paste0("tau_u_", var, "_", s_name_partition)
            tau_res_partition <- paste0("tau_res_", var)

            model_lines <- safe_add_lines(
                model_lines,
                c(
                    paste0("  # Variance Partitioning for independent var: ", var),
                    paste0("  ", sigma_total_param, " ~ dunif(0, 10)"),
                    paste0("  ", partition_param, " ~ dunif(0, 1)"),
                    paste0("  ", tau_total_param, " <- 1/(", sigma_total_param, " * ", sigma_total_param, ")"),
                    paste0("  ", tau_struct_partition, " <- ", tau_total_param, " / max(0.001, ", partition_param, ")"),
                    paste0("  ", tau_res_partition, " <- ", tau_total_param, " / max(0.001, 1 - ", partition_param, ")")
                )
            )
            # Map for post-processing
            param_map[[length(param_map) + 1]] <- list(response = var, predictor = s_name_partition, parameter = partition_param, type = "structure")
            param_map[[length(param_map) + 1]] <- list(response = var, predictor = s_name_partition, parameter = sigma_total_param, type = "structure")
            
        } else if (
          !is.null(fix_residual_variance) &&
            (var %in%
              names(fix_residual_variance) ||
              length(fix_residual_variance) == 1)
        ) {
          val <- if (var %in% names(fix_residual_variance)) {
            fix_residual_variance[[var]]
          } else {
            fix_residual_variance[[1]]
          }
          prec <- 1 / val # Inverse variance
          tau_line <- paste0(
            "  tau_res_",
            var,
            " <- ",
            prec,
            " # Fixed residual variance"
          )
          model_lines <- safe_add_lines(model_lines, tau_line)
        } else {
          tau_line <- paste0("  ", get_precision_prior(paste0("tau_res_", var), var))
          model_lines <- safe_add_lines(model_lines, tau_line)
        }
      } else {
        # [MODULAR PARTITIONING]
        # Check if partitioning was handled by a structure extension
        partition_already_handled <- FALSE
        s_name_handled <- NULL
        if (!is.null(structures)) {
          # The s_def in context is the last one from the loop (which would be the ONLY one if partitioning was used)
          if (exists("s_def") && isTRUE(s_def$partition_handled) && s_def$variable_name == var) {
            partition_already_handled <- TRUE
            # Heuristic: find which structure name was being used. 
            # Since partitioning only triggers if local_struct_count == 1, 
            # any structure sn that satisfies the mapping is the one.
            for (sn in names(structures)) {
                if (is_valid_structure_mapping(get_struct_lvl(sn, hierarchical_info), get_var_level(var, hierarchical_info), hierarchical_info, allow_identity = TRUE)) {
                    s_name_handled <- sn
                    break
                }
            }
          }
        }
        
        if (partition_already_handled) {
          # Register lambda in param_map
          param_map[[length(param_map) + 1]] <- list(
            response = var, 
            predictor = s_name_handled %||% "structure", 
            parameter = paste0("lambda_", var), 
            type = "structure"
          )
        } else {

          # FALLBACK: Multiple Structures or No Structure (Additive)
          if (
            !is.null(fix_residual_variance) &&
              (var %in%
                names(fix_residual_variance) ||
                length(fix_residual_variance) == 1)
          ) {
            val <- if (var %in% names(fix_residual_variance)) {
              fix_residual_variance[[var]]
            } else {
              fix_residual_variance[[1]]
            }
            prec <- 1 / val # Inverse variance
            tau_line <- paste0(
              "  tau_res_",
              var,
              " <- ",
              prec,
              " # Fixed residual variance"
            )
          } else {
            tau_line <- paste0("  ", get_precision_prior(paste0("tau_res_", var), var))
          }

          # Add residual variance if needed
          model_lines <- safe_add_lines(model_lines, tau_line)

          processed_ex_signals <- character(0)
          # [NOTE] Individual structures were already processed above in the loop (4148)
        }
      }

    }

    # Only generate TAU matrix with VCV when there is a structure defined
    if (FALSE && length(structure_names) > 0) {
      if (is_multi_structure) {
        model_lines <- c(
          model_lines,
          paste0(
            "  Mlam",
            var,
            " <- lambda",
            var,
            "*multiVCV[,,K] + (1 - lambda",
            var,
            ")*ID"
          ),
          paste0(
            "  TAU",
            tolower(var),
            " <- tau",
            var,
            "*inverse(Mlam",
            var,
            ")"
          )
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0(
            "  Mlam",
            var,
            " <- lambda",
            var,
            "*VCV + (1 - lambda",
            var,
            ")*ID"
          ),
          paste0(
            "  TAU",
            tolower(var),
            " <- tau",
            var,
            "*inverse(Mlam",
            var,
            ")"
          )
        )
      }
    }
  }

  # Verify if "ID" is actually used in the model (e.g. for residuals or GLMMs)
  # (Removed dummy usage of ID)

  # Ensure model is correctly terminated with a closing brace
  model_lines <- c(model_lines, "}")
  model_string <- paste(model_lines, collapse = "\n")

  # Convert param_map to data frame robustly
  param_map_df <- do.call(
    rbind,
    lapply(param_map, function(x) {
      if (!"type" %in% names(x)) x$type <- "coefficient"
      if (!"equation_index" %in% names(x)) x$equation_index <- NA
      # Ensure order and presence of all 5 fields
      fields <- c("response", "predictor", "parameter", "equation_index", "type")
      as.data.frame(x[fields], stringsAsFactors = FALSE)
    })
  )

  # [UNIFICATION] Master Registry Filter: Deduplicate Unified vs Legacy structural parameters
  # If both sigma_phylo_Var and sigma_Var_phylo exist, keep only sigma_phylo_Var
  if (!is.null(param_map_df) && nrow(param_map_df) > 0) {
      # Catch any sigma_Type_Variable parameters that might have legacy sigma_Variable_Type equivalents
      unified_names <- grep("^sigma_[a-zA-Z0-9]+_", param_map_df$parameter, value = TRUE)
      if (length(unified_names) > 0) {
          legacy_to_remove <- c()
          for (un in unified_names) {
              # Parse sigma_Structure_Variable -> Variable
              parts <- strsplit(un, "_")[[1]]
              if (length(parts) >= 3) {
                  s_type <- parts[2]
                  v_name <- paste(parts[3:length(parts)], collapse = "_")
                  legacy_pat <- paste0("sigma_", v_name, "_", s_type)
                  legacy_to_remove <- c(legacy_to_remove, legacy_pat)
              }
          }
          # Also remove corresponding tau legacy names
          legacy_tau_to_remove <- sub("sigma_", "tau_u_", legacy_to_remove)
          all_legacy <- c(legacy_to_remove, legacy_tau_to_remove)
          
          # Force Removal
          param_map_df <- param_map_df[!(param_map_df$parameter %in% all_legacy), ]
      }
  }

  return(list(model = model_string, parameter_map = param_map_df))
}
