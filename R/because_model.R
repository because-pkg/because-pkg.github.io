#' Generate a JAGS model string for Bayesian SEM (Because)
#'
#' This function builds the model code to be passed to JAGS based on a set of structural equations.
#' It supports custom covariance structures (spatial, phylogenetic, etc.).
#' Missing values are handled both in the response and predictor variables treating all of them as stochastic nodes.
#'
#' @param equations A list of model formulas.
#' @param is_multi_structure Logical; if \code{TRUE}, handles 3D Precision arrays (e.g. multi-object sampling).
#' @param variability Optional character vector or named character vector of variable names that have measurement error or within-species variability.
#' @param family Optional named character vector specifying the family/distribution for response variables.
#' @param vars_with_na Optional character vector of response variable names that have missing data.
#' @param induced_correlations Optional list of variable pairs with induced correlations from latent variables.
#' @param standardize_latent Logical (default TRUE). If TRUE, standardizes latent variables to unit variance.
#' @param structures A named list of structural objects (e.g. matrices, trees) to include as correlations.
#' @param latent Optional character vector of latent variable names.
#' @param poly_terms (Internal) List of polynomial terms for model generation.
#' @param fix_residual_variance Optional numeric value or named vector to fix residual variance.
#' @param latent_method Method for handling latent variables ("correlations" or "explicit").
#' @param random_structure_names Optional character vector of structural names applied to all variables.
#' @param random_terms Optional list of random effects (response, group).
#' @param categorical_vars Optional character vector of categorical variable names.
#' @param priors Optional named list of custom priors.
#' @param hierarchical_info Optional list containing multiscale data hierarchy (levels, link_vars).
#' @param engine Bayesian engine to use ("jags" or "nimble").
#' @param fix_latent Variable to fix for latent identification.
#' @param quiet Logical; suppress status messages.
#' @return A list with two elements:
#' \itemize{
#'   \item \code{model}: A character string containing the JAGS model code.
#'   \item \code{parameter_map}: A data frame mapping response variables to their predictors and parameter names.
#' }
#' @export
#' @importFrom stats formula terms setNames sd
#' @importFrom utils combn
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

  # Forward-compatibility: multi.tree used to be a standalone flag
  if (!is_multi_structure && "multi.tree" %in% names(match.call())) {
    is_multi_structure <- TRUE
  }

  # Helper: Get finest level for main loop N
  main_loop_N <- "N"
  if (!is.null(hierarchical_info)) {
    h_levels <- trimws(strsplit(hierarchical_info$hierarchy, ">")[[1]])
    finest_level <- h_levels[length(h_levels)]
    main_loop_N <- paste0("N_", finest_level)
  }

  # Populate random_structure_names from random_terms if not already provided
  if (length(random_terms) > 0) {
    rt_groups <- unique(sapply(random_terms, function(rt) rt$group))
    random_structure_names <- unique(c(random_structure_names, rt_groups))
  }

  # Deduplication: Purge random effects handled by structures
  if (!is.null(random_structure_names) && length(structure_names) > 0) {
    random_structure_names <- setdiff(random_structure_names, structure_names)
    
    if (!is.null(hierarchical_info)) {
      struct_lvls <- tolower(unique(na.omit(sapply(structure_names, function(s) get_struct_lvl(s, hierarchical_info)))))
      if (length(struct_lvls) > 0) {
        keep_r <- vapply(random_structure_names, function(r) {
          return(TRUE) 
        }, logical(1))
        random_structure_names <- random_structure_names[keep_r]
      }
    }
  }

  has_structure <- !is.null(structure_names) && length(structure_names) > 0
  has_random <- (!is.null(random_structure_names) && length(random_structure_names) > 0) ||
                (!is.null(random_terms) && length(random_terms) > 0)
  independent <- !has_structure && !has_random

  # Parse equations and extract deterministic terms
  deterministic_terms <- extract_deterministic_terms(equations)

  eq_list <- lapply(equations, function(eq) {
    response <- as.character(stats::formula(eq))[2]
    predictors <- attr(stats::terms(stats::formula(eq)), "term.labels")
    predictors <- sapply(predictors, sanitize_term_name)
    names(predictors) <- NULL
    list(response = response, predictors = predictors)
  })

  all_vars <- unique(unlist(lapply(eq_list, function(eq) {
    c(eq$response, eq$predictors)
  })))

  correlated_vars <- if (!is.null(induced_correlations)) {
    unique(unlist(induced_correlations))
  } else {
    character(0)
  }

  variability_list <- list()
  if (!is.null(variability)) {
    if (is.null(names(variability))) {
      variability_list <- stats::setNames(rep("se", length(variability)), variability)
    } else {
      variability_list <- variability
    }
  }

  dist_list <- list()
  if (!is.null(family)) {
    dist_list <- as.list(family)
  }

  # Build mutable context for compilation stages
  ctx <- new.env(parent = emptyenv())
  ctx$model_lines           <- character(0)
  ctx$param_map             <- list()
  ctx$declared_nodes        <- new.env(parent = emptyenv())
  ctx$declared_nodes$nodes  <- character(0)
  ctx$partitioned_responses <- character(0)
  ctx$vars_error_terms      <- list()
  ctx$beta_counter          <- list()
  ctx$response_counter      <- list()
  ctx$dist_list             <- dist_list
  ctx$eq_list               <- eq_list
  ctx$all_vars              <- all_vars
  ctx$correlated_vars       <- correlated_vars
  ctx$variability_list      <- variability_list
  ctx$main_loop_N           <- main_loop_N
  ctx$independent           <- independent
  ctx$structures            <- structures
  ctx$structure_names       <- structure_names
  ctx$random_structure_names<- random_structure_names
  ctx$random_terms          <- random_terms
  ctx$equations             <- equations
  ctx$latent                <- latent
  ctx$categorical_vars      <- categorical_vars
  ctx$priors                <- priors
  ctx$hierarchical_info     <- hierarchical_info
  ctx$family                <- family
  ctx$variability           <- variability
  ctx$vars_with_na          <- vars_with_na
  ctx$engine                <- engine
  ctx$quiet                 <- quiet
  ctx$is_multi_structure    <- is_multi_structure
  ctx$induced_correlations  <- induced_correlations
  ctx$fix_latent            <- fix_latent
  ctx$fix_residual_variance <- fix_residual_variance
  ctx$standardize_latent    <- standardize_latent
  ctx$poly_terms            <- poly_terms
  ctx$deterministic_terms   <- deterministic_terms

  # Execute modular pipeline stages sequentially
  generate_model_preamble(ctx)
  generate_linear_predictors(ctx)
  generate_likelihoods(ctx)
  generate_model_variability(ctx)
  generate_model_priors(ctx)

  # Terminate model with closing brace
  ctx$model_lines <- c(ctx$model_lines, "}")
  model_string <- paste(ctx$model_lines, collapse = "\n")

  # Convert param_map to data frame robustly
  param_map_df <- do.call(
    rbind,
    lapply(ctx$param_map, function(x) {
      if (!"type" %in% names(x)) x$type <- "coefficient"
      if (!"equation_index" %in% names(x)) x$equation_index <- NA
      fields <- c("response", "predictor", "parameter", "equation_index", "type")
      as.data.frame(x[fields], stringsAsFactors = FALSE)
    })
  )

  # Master Registry Filter: Deduplicate Unified vs Legacy structural parameters
  if (!is.null(param_map_df) && nrow(param_map_df) > 0) {
    unified_names <- grep("^sigma_[a-zA-Z0-9]+_", param_map_df$parameter, value = TRUE)
    if (length(unified_names) > 0) {
      legacy_to_remove <- c()
      for (un in unified_names) {
        parts <- strsplit(un, "_")[[1]]
        if (length(parts) >= 3) {
          s_type <- parts[2]
          v_name <- paste(parts[3:length(parts)], collapse = "_")
          legacy_pat <- paste0("sigma_", v_name, "_", s_type)
          legacy_to_remove <- c(legacy_to_remove, legacy_pat)
        }
      }
      legacy_tau_to_remove <- sub("sigma_", "tau_u_", legacy_to_remove)
      all_legacy <- c(legacy_to_remove, legacy_tau_to_remove)
      param_map_df <- param_map_df[!(param_map_df$parameter %in% all_legacy), ]
    }
  }

  return(list(model = model_string, parameter_map = param_map_df))
}
