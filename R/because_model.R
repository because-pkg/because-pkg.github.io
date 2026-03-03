#' Generate a JAGS model string for Phylogenetic Bayesian SEM (Because)
#'
#' This function builds the model code to be passed to JAGS based on a set of structural equations.
#' It supports both single and multiple phylogenetic trees (to account for phylogenetic uncertainty).
#' Missing values are handled both in the response and predictor variables treating all of them as stochastic nodes.
#'
#' @param equations A list of model formulas.
#' @param multi.tree Logical; if \code{TRUE}, incorporates phylogenetic uncertainty by sampling across a set of trees.
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
#'   \item Priors for intercepts (\code{alpha}), slopes (\code{beta}), lambda parameters (\code{lambda}), and residual precisions (\code{tau}).
#'   \item Phylogenetic covariance modeled via a single \code{VCV} matrix (when \code{multi.tree = FALSE}) or a 3D array \code{multiVCV[,,K]} with categorical sampling across trees (when \code{multi.tree = TRUE}).
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
#' cat(because_model(eqs, multi.tree = TRUE)$model)
#'
#' @param optimise Logical. If TRUE (default), use random effects formulation for 4.6× speedup.
#'   If FALSE, use original marginal covariance formulation.
#' @param standardize_latent Logical (default TRUE). If TRUE, standardizes latent variables to unit variance.
#' @param structure_names (Internal) Character vector of names for multiple trees/structures.
#' @param latent Optional character vector of latent variable names.
#' @param poly_terms (Internal) List of polynomial terms for model generation.
#' @param fix_residual_variance Optional numeric value or named vector to fix residual variance.
#' @export
#' @importFrom stats formula terms setNames sd
#' @importFrom utils combn
#'
because_model <- function(
  equations,
  multi.tree = FALSE,
  latent_method = "correlations",
  structure_names = "phylo",
  structures = NULL,
  random_structure_names = NULL,
  random_terms = list(),
  vars_with_na = NULL,
  induced_correlations = NULL,
  variability = NULL,
  family = NULL,
  optimise = TRUE,
  standardize_latent = TRUE,
  poly_terms = NULL,
  latent = NULL,
  categorical_vars = NULL,
  fix_residual_variance = NULL,
  priors = NULL,
  hierarchical_info = NULL
) {
  # Helper: returns b if a is NULL or if a is a list element that doesn't exist
  `%||%` <- function(a, b) {
    tryCatch(if (!is.null(a)) a else b, error = function(e) b)
  }

  # Helper: Get prior for a parameter (custom override or default)
  get_prior <- function(param_name, default_prior) {
    if (!is.null(priors) && param_name %in% names(priors)) {
      return(paste0(param_name, " ~ ", priors[[param_name]]))
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
    for (lvl_name in names(h_info$levels)) {
      if (var %in% h_info$levels[[lvl_name]]) {
        return(lvl_name)
      }
    }

    return(NULL) # Variable not found in any level
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

  # Helper: Get level for a structure
  get_struct_lvl <- function(s_name, h_info) {
    if (is.null(h_info) || is.null(h_info$structure_levels)) {
      return(NULL)
    }
    return(h_info$structure_levels[[s_name]])
  }

  # Helper: Get index for mapping structure levels to response levels
  get_struct_index <- function(s_name, response, h_info) {
    s_lvl <- get_struct_lvl(s_name, h_info)
    r_lvl <- get_var_level(response, h_info)
    if (is.null(s_lvl) || is.null(r_lvl) || s_lvl == r_lvl) {
      return("i")
    }
    # Level mismatch: use the indicator mapping e.g. site_idx_obs[i]
    return(paste0(s_lvl, "_idx_", r_lvl, "[i]"))
  }

  # Helper: Validate if a structure's level can map to a response's level
  is_valid_structure_mapping <- function(s_lvl, r_lvl, h_info) {
    if (is.null(s_lvl) || is.null(r_lvl) || s_lvl == r_lvl) {
      return(TRUE)
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

  # Helper: Get index expression for accessing a predictor from a coarser level
  # Returns "pred[Coarser_idx_Finer[i]]" or "pred[i]" if same level
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

  # Helper: Check if random effect grouping level is compatible with response level
  # Valid if: Rank(Group) <= Rank(Response) (Coarser or Equal to Response)
  is_valid_random_level <- function(response, group_var, h_info) {
    if (is.null(h_info)) {
      return(TRUE)
    }

    resp_lvl <- get_var_level(response, h_info)
    grp_lvl <- get_var_level(group_var, h_info)

    if (is.null(resp_lvl) || is.null(grp_lvl)) {
      return(TRUE)
    } # Assume true if unknown

    # Use the robust checker
    return(is_valid_structure_mapping(grp_lvl, resp_lvl, h_info))
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
        if (!is_valid_structure_mapping(v_lvl, lvl, h_info)) {
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

    # Return the coarsest of the descendants (the intersection)
    coarsest <- candidates[1]
    if (length(candidates) > 1) {
      for (i in 2:length(candidates)) {
        if (is_valid_structure_mapping(candidates[i], coarsest, h_info)) {
          coarsest <- candidates[i]
        }
      }
    }
    return(coarsest)
  }

  # Compute main loop N once for use throughout
  main_loop_N <- "N"
  if (!is.null(hierarchical_info)) {
    h_levels <- trimws(strsplit(hierarchical_info$hierarchy, ">")[[1]])
    finest_level <- h_levels[length(h_levels)]
    main_loop_N <- paste0("N_", finest_level)
  }

  # Helper: Get N string (for hierarchical-aware array dimensions)
  N_str <- function() main_loop_N

  has_structure <- !is.null(structure_names) && length(structure_names) > 0
  has_random <- !is.null(random_structure_names) &&
    length(random_structure_names) > 0
  independent <- !has_structure && !has_random

  # Flag for multi-structure (e.g., multiPhylo with 3D precision arrays)
  # This applies to ANY structure type when multi.tree is TRUE
  is_multi_structure <- multi.tree

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
  # Check if we need zero_vec (for induced_correlations OR multinomial)
  need_zero_vec <- !is.null(induced_correlations)
  if (!need_zero_vec && !is.null(family)) {
    if (any(family == "multinomial") || any(family == "ordinal")) {
      need_zero_vec <- TRUE
    }
  }

  model_lines <- c(
    "model {",
    "  # Common structures and priors"
  )

  if (need_zero_vec) {
    model_lines <- c(
      model_lines,
      "  for(k in 1:N) { zero_vec[k] <- 0 }",
      "  zero_vec_2[1] <- 0",
      "  zero_vec_2[2] <- 0"
    )
  }

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

  model_lines <- c(model_lines, "  # Structural equations")

  # --- Generic Structure Setup ---
  if (!is.null(structures)) {
    for (s_name in names(structures)) {
      # Dispatch to S3 generic
      def <- jags_structure_definition(
        structures[[s_name]],
        optimize = optimise
      )
      if (!is.null(def$setup_code)) {
        model_lines <- c(model_lines, def$setup_code)
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

    # Legacy Multi-tree Setup (Fallback if structures not provided)
    if (multi.tree && is.null(structures)) {
      model_lines <- c(
        model_lines,
        "  # Multi-tree sampling",
        "  K ~ dcat(p_tree[])",
        "  for (k in 1:Ntree) {",
        "    p_tree[k] <- 1/Ntree",
        "  }"
      )
    }

    if (length(exogenous_vars) > 0) {
      model_lines <- c(
        model_lines,
        "  # Priors for exogenous latent variables (variable with error but no parent)",
        "  for (i in 1:N) {"
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
        vars_in_expr <- all.vars(parse(text = dt$original))
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
  model_lines <- c(model_lines, "")

  for (j in seq_along(eq_list)) {
    eq <- eq_list[[j]]
    response <- eq$response
    predictors <- eq$predictors
    dist <- dist_list[[response]] %||% "gaussian"

    # Hierarchical Loop Wrapper
    eq_loop_N <- get_loop_bound(response, hierarchical_info)
    model_lines <- c(model_lines, paste0("  for (i in 1:", eq_loop_N, ") {"))

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
      model_lines <- c(model_lines, paste0("    ", response, "[i] <- ", expr))
      # No parameters to monitor or priors to add for this equation
      next
    }

    # Count and assign unique suffix for the response variable
    response_count <- response_counter[[response]] %||% 0
    response_count <- response_count + 1
    response_counter[[response]] <- response_count
    suffix <- if (response_count == 1) "" else as.character(response_count)

    alpha <- paste0("alpha_", response, suffix)
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
        equation_index = j
      )
    }

    if (dist == "gaussian" || dist == "occupancy" || grepl("^p_", response)) {
      mu <- paste0("mu_", response, suffix)
      model_lines <- c(model_lines, paste0("    ", mu, "[i] <- ", linpred))
    } else if (dist == "binomial") {
      # Binomial: logit(p) = linpred + error
      # error ~ dmnorm(0, TAU)
      # We define the error mean here as 0
      mu_err <- paste0("mu_err_", response, suffix)
      err <- paste0("err_", response, suffix)
      p <- paste0("p_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    ", mu_err, "[i] <- 0"),
        paste0("    logit(", p, "[i]) <- ", linpred, " + ", err, "[i]"),
        paste0("    ", response, "[i] ~ dbern(", p, "[i])"),
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
        linpred_k <- paste0(linpred_k, " + ", beta_name, "[k] * ", pred, "[i]")

        # Map (only once)
        key <- paste(response, pred, suffix, sep = "_")
        if (!key %in% names(beta_counter)) {
          beta_counter[[key]] <- beta_name # Mark as used
          param_map[[length(param_map) + 1]] <- list(
            response = response,
            predictor = pred,
            parameter = paste0(beta_name, "[]"),
            equation_index = j
          )
        }
      }

      linpred_k <- paste0(linpred_k, " + ", err, "[i, k]")

      model_lines <- c(
        model_lines,
        paste0("      L_", response, "[i, k] <- ", linpred_k),
        "    }"
      )

      # Softmax
      model_lines <- c(
        model_lines,
        paste0("    # Softmax for ", response),
        paste0("    for (k in 1:", K_var, ") {"),
        paste0(
          "      exp_L_",
          response,
          "[i, k] <- exp(L_",
          response,
          "[i, k])"
        ),
        "    }",
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
          "[i, k] <- exp_L_",
          response,
          "[i, k] / sum_exp_L_",
          response,
          "[i]"
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
        ),
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
    } else if (dist == "ordinal") {
      # Ordinal: Cumulative Logit (Proportional Odds)
      # P(Y <= k) = logit^(-1)(cutpoint[k] - eta)
      # eta = linpred + error

      K_var <- paste0("K_", response)
      err <- paste0("err_", response, suffix)
      eta <- paste0("eta_", response, suffix)

      # Linear predictor (eta)
      # Note: No intercept in eta (intercept is absorbed into cutpoints)
      linpred_no_int <- "0"
      for (pred in predictors) {
        beta_name <- paste0("beta_", response, "_", pred)
        linpred_no_int <- paste0(
          linpred_no_int,
          " + ",
          beta_name,
          " * ",
          pred,
          "[i]"
        )

        # Map
        key <- paste(response, pred, suffix, sep = "_")
        if (!key %in% names(beta_counter)) {
          beta_counter[[key]] <- beta_name
          param_map[[length(param_map) + 1]] <- list(
            response = response,
            predictor = pred,
            parameter = beta_name,
            equation_index = j
          )
        }
      }

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
        ),
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
    } else if (dist == "poisson") {
      # Poisson: log(μ) = linpred + error
      # Naturally handles overdispersion via epsilon

      err <- paste0("err_", response, suffix)
      mu <- paste0("mu_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    # Poisson log link for ", response),
        paste0("    log(", mu, "[i]) <- ", linpred, " + ", err, "[i]"),
        paste0("    ", response, "[i] ~ dpois(", mu, "[i])"),
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
    } else if (dist == "negbinomial") {
      # Negative Binomial: log(μ) = linpred + error
      # Y ~ NegBin(p, r) where p = r/(r+μ) and r = size parameter

      err <- paste0("err_", response, suffix)
      mu <- paste0("mu_", response, suffix)
      p <- paste0("p_", response, suffix)
      r <- paste0("r_", response, suffix)

      model_lines <- c(
        model_lines,
        paste0("    # Negative Binomial log link for ", response),
        paste0("    log(", mu, "[i]) <- ", linpred, " + ", err, "[i]"),
        paste0("    ", p, "[i] <- ", r, " / (", r, " + ", mu, "[i])"),
        paste0("    ", response, "[i] ~ dnegbin(", p, "[i], ", r, ")"),
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
        paste0("    log(", mu, "[i]) <- ", linpred, " + ", err, "[i]"),

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
          ") * exp(logdensity.pois(",
          response,
          "[i], ",
          mu,
          "[i]))"
        ),

        # Select likelihood based on Y[i]
        # Use a small epsilon to avoid log(0)
        paste0(
          "    lik_",
          response,
          "[i] <- ifelse(",
          response,
          "[i] == 0, lik_zero_",
          response,
          "[i], lik_pos_",
          response,
          "[i])"
        ),

        paste0(
          "    log_lik_",
          response,
          suffix,
          "[i] <- log(lik_",
          response,
          "[i])"
        ),
        paste0(
          "    phi_",
          response,
          "[i] <- -log_lik_",
          response,
          suffix,
          "[i] + 10000"
        ),
        paste0("    zeros[i] ~ dpois(phi_", response, "[i])")
      )
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
        paste0("    log(", mu, "[i]) <- ", linpred, " + ", err, "[i]"),
        paste0("    ", p, "[i] <- ", r, " / (", r, " + ", mu, "[i])"),

        # Zeros trick
        # Zero case: psi + (1-psi) * p^r
        paste0(
          "    lik_zero_",
          response,
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
          "[i] <- (1-",
          psi,
          ") * exp(logdensity.negbin(",
          response,
          "[i], ",
          p,
          "[i], ",
          r,
          "))"
        ),

        # Select likelihood
        paste0(
          "    lik_",
          response,
          "[i] <- ifelse(",
          response,
          "[i] == 0, lik_zero_",
          response,
          "[i], lik_pos_",
          response,
          "[i])"
        ),

        paste0(
          "    log_lik_zero_",
          response,
          "[i] <- log(max(1.0E-30, ",
          "lik_zero_",
          response,
          "[i]))"
        ),
        paste0(
          "    log_lik_pos_",
          response,
          "[i] <- log(max(1.0E-30, 1-",
          psi,
          ")) + ",
          "logdensity.negbin(",
          response,
          "[i], ",
          p,
          "[i], ",
          r,
          ")"
        ),
        paste0(
          "    log_lik_",
          response,
          suffix,
          "[i] <- ifelse(",
          response,
          "[i] == 0, log_lik_zero_",
          response,
          "[i], log_lik_pos_",
          response,
          "[i])"
        ),
        paste0(
          "    phi_",
          response,
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
    } else if (dist == "occupancy") {
      # No Zeros-trick or custom likelihood loop needed here
      # The occupancy likelihood is handled in the Likelihoods section.
    } else {
      stop(paste("Unknown distribution:", dist))
    }

    # Close Hierarchical Loop
    model_lines <- c(model_lines, "  }")
  }

  model_lines <- c(model_lines, "  # Multivariate normal likelihoods")

  # Likelihoods for responses
  for (response in names(response_counter)) {
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
        if (!is.null(vars_with_na) && response %in% vars_with_na && !optimise) {
          # Use Latent Variable (GLMM) approach for missing data (Only if optimisation disabled)
          # Y[i] ~ dnorm(mu[i] + err[i], tau_res)
          # err[1:N] ~ dmnorm(0, tau_phylo * inv(VCV))

          err <- paste0("err_", response, suffix)
          mu_err <- paste0("mu_err_", response, suffix)
          tau_res <- paste0("tau_res_", response, suffix)

          model_lines <- c(
            model_lines,
            paste0(
              "  # GLMM likelihood for missing data (preserves phylo signal)"
            ),
            paste0("  for (i in 1:", loop_bound, ") {"),
            paste0(
              "    ",
              response_var,
              "[i] ~ dnorm(",
              mu,
              "[i] + ",
              err,
              "[i], ",
              tau_res,
              ")"
            )
          )
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
              "[i] + ",
              err,
              "[i], ",
              tau_res,
              ")"
            ),
            paste0("  }")
          )
        } else {
          if (independent) {
            # Independent Model (No random effects)
            # Y[i] ~ dnorm(mu[i], tau)
            # Note: We use tau_e for consistency with optimized model naming if desired,
            # but tau is standard for simple normal. Let's use tau_e to distinguish from matrix TAU
            tau_e <- paste0("tau_e_", response, suffix)

            model_lines <- c(
              model_lines,
              paste0("  for (i in 1:", loop_bound, ") {"),
              paste0(
                "    ",
                response_var,
                "[i] ~ dnorm(",
                mu,
                "[i], ",
                tau_e,
                ")"
              )
            )
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
                tau_e,
                ")"
              ),
              paste0("  }")
            )
          } else if (optimise) {
            # Optimized Random Effects Formulation (Additive)
            additive_terms <- ""

            for (s_name in structure_names) {
              if (
                !is_valid_structure_mapping(
                  get_struct_lvl(s_name, hierarchical_info),
                  get_var_level(response, hierarchical_info),
                  hierarchical_info
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
              u_std <- paste0("u_std_", response, suffix, s_suffix)
              u <- paste0("u_", response, suffix, s_suffix)
              tau_u <- paste0("tau_u_", response, suffix, s_suffix)

              s_lvl <- get_struct_lvl(s_name, hierarchical_info)
              s_bound <- if (is.null(s_lvl)) loop_bound else paste0("N_", s_lvl)
              s_zeros <- if (is.null(s_lvl)) {
                "zeros"
              } else {
                paste0("zeros_", s_lvl)
              }

              # Call Generic to define the error term u ~ dmnorm(...)
              def <- jags_structure_definition(
                structures[[s_name]],
                variable_name = u_std,
                optimize = optimise,
                precision_parameter = tau_u
              )

              if (!is.null(def$setup_code)) {
                model_lines <- c(model_lines, def$setup_code)
              }

              prec_name <- paste0("Prec_", s_name)
              prec_index <- if (!is.null(def$prec_index)) {
                def$prec_index
              } else if (multi.tree && is_multi_structure) {
                paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, ", K]")
              } else {
                paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, "]")
              }

              model_lines <- c(
                model_lines,
                paste0(
                  "  ",
                  u_std,
                  "[1:",
                  s_bound,
                  "] ~ dmnorm(",
                  s_zeros,
                  "[1:",
                  s_bound,
                  "], ",
                  prec_index,
                  ")"
                ),
                paste0(
                  "  for (i in 1:",
                  loop_bound,
                  ") { ",
                  u,
                  "[i] <- ",
                  u_std,
                  "[",
                  get_struct_index(s_name, response, hierarchical_info),
                  "] / sqrt(",
                  tau_u,
                  ") }"
                )
              )
              additive_terms <- paste0(additive_terms, " + ", u, "[i]")
            }

            # Random Effects (Grouped)
            for (r_name in random_structure_names) {
              # Check relevance using random_terms mapping
              is_relevant <- TRUE
              if (length(random_terms) > 0) {
                matches <- Filter(
                  function(x) x$response == response && x$group == r_name,
                  random_terms
                )
                if (length(matches) == 0) is_relevant <- FALSE
              }

              if (is_relevant) {
                # Valid hierarchical check: Group <= Response (Coarser or Equal)
                if (
                  !is_valid_random_level(response, r_name, hierarchical_info)
                ) {
                  is_relevant <- FALSE
                }
              }

              if (is_relevant) {
                s_suffix <- paste0("_", r_name)

                u_std <- paste0("u_std_", response, suffix, s_suffix)
                u <- paste0("u_", response, suffix, s_suffix)
                tau_u <- paste0("tau_u_", response, suffix, s_suffix)

                n_groups <- get_loop_bound(r_name, hierarchical_info)

                # Smart Group Index Logic for Hierarchy
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

                  # If levels differ, use nested transitive index
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
                prec_name <- paste0("Prec_", r_name)
                zeros_name <- paste0("zeros_", r_name)

                model_lines <- c(
                  model_lines,
                  paste0(
                    "  ",
                    u_std,
                    "[1:",
                    n_groups,
                    "] ~ dmnorm(",
                    zeros_name,
                    "[1:",
                    n_groups,
                    "], ",
                    prec_name,
                    "[1:",
                    n_groups,
                    ", 1:",
                    n_groups,
                    "])"
                  ),
                  paste0("  for (g in 1:", n_groups, ") {"),
                  paste0(
                    "    ",
                    u,
                    "[g] <- ",
                    u_std,
                    "[g] / sqrt(",
                    tau_u,
                    ")"
                  ),
                  paste0("  }")
                )

                additive_terms <- paste0(
                  additive_terms,
                  " + ",
                  u,
                  "[",
                  group_idx,
                  "]"
                )
              }
            } # End relevant check

            tau_e <- paste0("tau_e_", response, suffix)

            model_lines <- c(
              model_lines,
              paste0("  for (i in 1:", loop_bound, ") {"),
              paste0(
                "    ",
                response_var,
                "[i] ~ dnorm(",
                mu,
                "[i]",
                additive_terms,
                ", ",
                tau_e,
                ")"
              )
            )
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
                tau_e,
                ")"
              ),
              paste0("  }")
            )
          } else {
            # Standard MVN for complete data (Marginal)
            # Note: dmnorm is joint, so we compute pointwise log-lik separately
            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                response_var,
                "[1:",
                loop_bound,
                "] ~ dmnorm(",
                mu,
                "[1:",
                loop_bound,
                "], ",
                tau,
                ")"
              )
            )
            model_lines <- c(
              model_lines,
              "  # Pointwise log-likelihood for MVN",
              paste0("  for (i in 1:", loop_bound, ") {"),
              paste0(
                "    tau_marg_",
                response,
                suffix,
                "[i] <- ",
                tau,
                "[i, i]  # Extract diagonal precision (marginal variance)"
              ),
              paste0(
                "    log_lik_",
                response,
                suffix,
                "[i] <- logdensity.norm(",
                response_var,
                "[i], ",
                mu,
                "[i], tau_marg_",
                response,
                suffix,
                "[i])"
              ),
              paste0("  }")
            )
          }
        }
      } else if (dist == "binomial") {
        # For binomial, the error term has the phylogenetic structure
        err <- paste0("err_", response, suffix)

        if (independent) {
          # Independent Binomial: Standard GLM (no residual error)
          model_lines <- c(
            model_lines,
            paste0("  # Independent (Standard GLM) for binomial: ", response),
            paste0("  for (i in 1:N) {"),
            paste0("    ", err, "[i] <- 0"),
            paste0("  }")
          )
        } else if (optimise) {
          # Optimized Random Effects Formulation
          epsilon <- paste0("epsilon_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          # Initialize accumulator for random effects
          total_u <- ""

          # Phylogenetic / N-dim Structures
          for (s_name in structure_names) {
            if (
              !is_valid_structure_mapping(
                get_struct_lvl(s_name, hierarchical_info),
                get_var_level(response, hierarchical_info),
                hierarchical_info
              )
            ) {
              next
            }
            # Structure properties
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_suffix <- paste0("_", s_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            # Call Generic
            def <- jags_structure_definition(
              structures[[s_name]],
              variable_name = u_std,
              optimize = optimise,
              precision_parameter = tau_u
            )

            if (!is.null(def$setup_code)) {
              model_lines <- c(model_lines, def$setup_code)
            }

            prec_name <- paste0("Prec_", s_name)
            prec_index <- if (!is.null(def$prec_index)) {
              def$prec_index
            } else if (multi.tree && is_multi_structure) {
              paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, ", K]")
            } else {
              paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, "]")
            }

            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                u_std,
                "[1:",
                s_bound,
                "] ~ dmnorm(",
                s_zeros,
                "[1:",
                s_bound,
                "], ",
                prec_index,
                ")"
              ),
              paste0(
                "  for (i in 1:N) { ",
                u,
                "[i] <- ",
                u_std,
                "[",
                get_struct_index(s_name, response, hierarchical_info),
                "] / sqrt(",
                tau_u,
                ") }"
              )
            )
            total_u <- paste0(total_u, " + ", u, "[i]")
          }

          # Random Group Structures
          for (r_name in random_structure_names) {
            # Valid hierarchical check: Group <= Response
            if (!is_valid_random_level(response, r_name, hierarchical_info)) {
              next
            }

            # Calculate loop bound for response (needed for Binomial/GLMM loop)
            loop_bound <- get_loop_bound(response, hierarchical_info)
            s_suffix <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            n_groups <- paste0("N_", r_name)

            # Smart Group Index Logic for Hierarchy
            group_idx <- paste0("group_", r_name, "[i]")
            if (!is.null(hierarchical_info)) {
              resp_bound <- get_loop_bound(response, hierarchical_info)
              grp_bound <- get_loop_bound(r_name, hierarchical_info)

              # Strip N_ prefix
              resp_lvl <- sub("^N_", "", resp_bound)
              grp_lvl <- sub("^N_", "", grp_bound)

              # Resolve "N" to finest level name
              finest_lvl <- trimws(strsplit(hierarchical_info$hierarchy, ">")[[
                1
              ]])
              finest_lvl <- finest_lvl[length(finest_lvl)]

              if (resp_lvl == "N") {
                resp_lvl <- finest_lvl
              }
              if (grp_lvl == "N") {
                grp_lvl <- finest_lvl
              }

              # If levels differ, use nested transitive index
              # This creates u[ group_var[ bridge_idx[i] ] ]
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
            prec_name <- paste0("Prec_", r_name)
            zeros_name <- paste0("zeros_", r_name)

            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                u_std,
                "[1:",
                n_groups,
                "] ~ dmnorm(",
                zeros_name,
                "[1:",
                n_groups,
                "], ",
                prec_name,
                "[1:",
                n_groups,
                ", 1:",
                n_groups,
                "])"
              ),
              paste0("  for (g in 1:", n_groups, ") {"),
              paste0("    ", u, "[g] <- ", u_std, "[g] / sqrt(", tau_u, ")"),
              paste0("  }")
            )
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
          }

          model_lines <- c(
            model_lines,
            paste0("  # Binomial error term summation: ", response),
            paste0("  for (i in 1:", loop_bound, ") {"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")")
          )

          # Sum all terms: overdispersion + random effects + MAG terms
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) {
            paste0(" + ", paste(mag_terms, collapse = " + "))
          } else {
            ""
          }

          model_lines <- c(
            model_lines,
            paste0("    ", err, "[i] <- ", epsilon, "[i]", total_u, total_mag),
            paste0("  }")
          )
        } else {
          # Marginal Formulation (original)
          mu_err <- paste0("mu_err_", response, suffix)
          model_lines <- c(
            model_lines,
            paste0("  ", err, "[1:N] ~ dmnorm(", mu_err, "[], ", tau, ")")
          )
        }
      } else if (dist == "multinomial") {
        # Multinomial error terms: err[1:N, k]
        # Independent phylogenetic effects for each k (2..K)
        err <- paste0("err_", response, suffix)
        K_var <- paste0("K_", response)

        if (independent) {
          # Independent Multinomial: Standard GLM
          model_lines <- c(
            model_lines,
            paste0(
              "  # Independent (Standard GLM) for multinomial: ",
              response
            ),
            paste0("  for (k in 2:", K_var, ") {"),
            paste0("    for (i in 1:N) {"),
            paste0("      ", err, "[i, k] <- 0"),
            paste0("    }"),
            paste0("  }")
          )
        } else if (optimise) {
          # Optimised Random Effects Formulation for Multinomial
          epsilon <- paste0("epsilon_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for multinomial: ", response),
            paste0("  for (k in 2:", K_var, ") {")
          )

          # Initialize accumulator for this category k
          total_u <- ""

          # Phylogenetic / N-dim Structures
          for (s_name in structure_names) {
            if (
              !is_valid_structure_mapping(
                get_struct_lvl(s_name, hierarchical_info),
                get_var_level(response, hierarchical_info),
                hierarchical_info
              )
            ) {
              next
            }
            # Structure properties
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_suffix <- paste0("_", s_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            prec_name <- paste0("Prec_", s_name)
            prec_index <- if (multi.tree && is_multi_structure) {
              paste0(prec_name, "[1:N, 1:N, k]")
            } else {
              paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, "]")
            }

            model_lines <- c(
              model_lines,
              paste0(
                "    ",
                u_std,
                "[1:",
                s_bound,
                ", k] ~ dmnorm(",
                s_zeros,
                "[1:",
                s_bound,
                "], ",
                prec_index,
                ")"
              ),
              paste0(
                "    for (i in 1:",
                get_loop_bound(response, hierarchical_info),
                ") { ",
                u,
                "[i, k] <- ",
                u_std,
                "[",
                get_struct_index(s_name, response, hierarchical_info),
                ", k] / sqrt(",
                tau_u,
                "[k]) }"
              )
            )
            total_u <- paste0(total_u, " + ", u, "[i, k]")
          }

          # Random Group Structures
          for (r_name in random_structure_names) {
            s_suffix <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            n_groups <- paste0("N_", r_name)

            # Smart Group Index Logic for Hierarchy
            group_idx <- paste0("group_", r_name, "[i]")
            if (!is.null(hierarchical_info)) {
              resp_bound <- get_loop_bound(response, hierarchical_info)
              grp_bound <- get_loop_bound(r_name, hierarchical_info)

              # Strip N_ prefix
              resp_lvl <- sub("^N_", "", resp_bound)
              grp_lvl <- sub("^N_", "", grp_bound)

              # Resolve "N" to finest level name
              finest_lvl <- trimws(strsplit(hierarchical_info$hierarchy, ">")[[
                1
              ]])
              finest_lvl <- finest_lvl[length(finest_lvl)]

              if (resp_lvl == "N") {
                resp_lvl <- finest_lvl
              }
              if (grp_lvl == "N") {
                grp_lvl <- finest_lvl
              }

              # If levels differ, use nested transitive index
              # This creates u[ group_var[ bridge_idx[i] ] ]
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
            prec_name <- paste0("Prec_", r_name)
            zeros_name <- paste0("zeros_", r_name)

            model_lines <- c(
              model_lines,
              paste0(
                "    ",
                u_std,
                "[1:",
                n_groups,
                ", k] ~ dmnorm(",
                zeros_name,
                "[1:",
                n_groups,
                "], ",
                prec_name,
                "[1:",
                n_groups,
                ", 1:",
                n_groups,
                "])"
              ),
              paste0("    for (g in 1:", n_groups, ") {"),
              paste0(
                "      ",
                u,
                "[g, k] <- ",
                u_std,
                "[g, k] / sqrt(",
                tau_u,
                "[k])"
              ),
              paste0("    }")
            )

            total_u <- paste0(total_u, " + ", u, "[", group_idx, ", k]")
          }

          model_lines <- c(
            model_lines,
            paste0("    for (i in 1:N) {"),
            paste0("      ", epsilon, "[i, k] ~ dnorm(0, ", tau_e, "[k])"),
            paste0("      ", err, "[i, k] <- ", epsilon, "[i, k]", total_u),
            paste0("    }"),
            paste0("  }")
          )
        } else {
          model_lines <- c(
            model_lines,
            paste0("  # Multinomial phylogenetic errors for ", response),
            paste0("  for (k in 2:", K_var, ") {"),
            paste0(
              "    ",
              err,
              "[1:N, k] ~ dmnorm(zero_vec[], TAU_",
              tolower(response),
              "_",
              suffix,
              "[,,k])"
            ),
            "  }"
          )
        }
      } else if (dist == "ordinal") {
        # Ordinal error term: err[1:N]
        # Single phylogenetic effect (unlike multinomial with K-1 effects)
        err <- paste0("err_", response, suffix)

        if (independent) {
          # Independent Ordinal: Standard GLM
          model_lines <- c(
            model_lines,
            paste0("  # Independent (Standard GLM) for ordinal: ", response),
            paste0("  for (i in 1:N) {"),
            paste0("    ", err, "[i] <- 0"),
            paste0("  }")
          )
        } else {
          # Random Effects Formulation (Default/Optimised)
          # Note: Ordinal currently only supports optimized formulation in this codebase structure
          u_std <- paste0("u_std_", response, suffix)
          u <- paste0("u_", response, suffix)
          epsilon <- paste0("epsilon_", response, suffix)
          tau_u <- paste0("tau_u_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          # Handle multi-tree with generic structure naming
          s_name <- if (length(structure_names) > 0) {
            structure_names[1]
          } else {
            "struct"
          }

          s_lvl <- get_struct_lvl(s_name, hierarchical_info)
          s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
          s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)

          prec_name <- paste0("Prec_", s_name)
          prec_index <- if (is_multi_structure) {
            paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, ", K]")
          } else {
            paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, "]")
          }

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for ordinal: ", response),
            paste0(
              "  ",
              u_std,
              "[1:",
              s_bound,
              "] ~ dmnorm(",
              s_zeros,
              "[1:",
              s_bound,
              "], ",
              prec_index,
              ")"
            )
          )
          # Sum MAG terms
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) {
            paste0(" + ", paste(mag_terms, collapse = " + "))
          } else {
            ""
          }

          model_lines <- c(
            model_lines,
            paste0("  for (i in 1:N) {"),
            paste0(
              "    ",
              u,
              "[i] <- ",
              u_std,
              "[",
              get_struct_index(s_name, response, hierarchical_info),
              "] / sqrt(",
              tau_u,
              ")"
            ),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0(
              "    ",
              err,
              "[i] <- ",
              u,
              "[i] + ",
              epsilon,
              "[i]",
              total_mag
            ),
            paste0("  }")
          )
        }
      } else if (dist == "poisson" || dist == "zip") {
        # Poisson error term: err[1:N]
        # Single phylogenetic effect (like ordinal)
        err <- paste0("err_", response, suffix)

        if (independent) {
          # Independent Poisson: Standard GLM (no overdispersion unless MAG present)
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) {
            paste0(" + ", paste(mag_terms, collapse = " + "))
          } else {
            ""
          }

          model_lines <- c(
            model_lines,
            paste0("  # Independent (Standard GLM) for Poisson: ", response),
            paste0("  for (i in 1:N) {"),
            paste0("    ", err, "[i] <- 0", total_mag),
            paste0("  }")
          )
        } else if (optimise) {
          # Optimized Random Effects
          epsilon <- paste0("epsilon_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          # Initialize accumulator
          total_u <- ""

          # Phylogenetic / N-dim Structures
          for (s_name in structure_names) {
            if (
              !is_valid_structure_mapping(
                get_struct_lvl(s_name, hierarchical_info),
                get_var_level(response, hierarchical_info),
                hierarchical_info
              )
            ) {
              next
            }
            # Structure properties
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_suffix <- paste0("_", s_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            prec_name <- paste0("Prec_", s_name)
            prec_index <- if (multi.tree && is_multi_structure) {
              paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, ", K]")
            } else {
              paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, "]")
            }

            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                u_std,
                "[1:",
                s_bound,
                "] ~ dmnorm(",
                s_zeros,
                "[1:",
                s_bound,
                "], ",
                prec_index,
                ")"
              ),
              paste0(
                "  for (i in 1:",
                get_loop_bound(response, hierarchical_info),
                ") { ",
                u,
                "[i] <- ",
                u_std,
                "[",
                get_struct_index(s_name, response, hierarchical_info),
                "] / sqrt(",
                tau_u,
                ") }"
              )
            )
            total_u <- paste0(total_u, " + ", u, "[i]")
          }

          # Random Group Structures
          for (r_name in random_structure_names) {
            s_suffix <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            n_groups <- paste0("N_", r_name)

            # Smart Group Index Logic for Hierarchy
            group_idx <- paste0("group_", r_name, "[i]")
            if (!is.null(hierarchical_info)) {
              resp_bound <- get_loop_bound(response, hierarchical_info)
              grp_bound <- get_loop_bound(r_name, hierarchical_info)

              # Strip N_ prefix
              resp_lvl <- sub("^N_", "", resp_bound)
              grp_lvl <- sub("^N_", "", grp_bound)

              # Resolve "N" to finest level name
              finest_lvl <- trimws(strsplit(hierarchical_info$hierarchy, ">")[[
                1
              ]])
              finest_lvl <- finest_lvl[length(finest_lvl)]

              if (resp_lvl == "N") {
                resp_lvl <- finest_lvl
              }
              if (grp_lvl == "N") {
                grp_lvl <- finest_lvl
              }

              # If levels differ, use nested transitive index
              # This creates u[ group_var[ bridge_idx[i] ] ]
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
            prec_name <- paste0("Prec_", r_name)
            zeros_name <- paste0("zeros_", r_name)

            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                u_std,
                "[1:",
                n_groups,
                "] ~ dmnorm(",
                zeros_name,
                "[1:",
                n_groups,
                "], ",
                prec_name,
                "[1:",
                n_groups,
                ", 1:",
                n_groups,
                "])"
              ),
              paste0("  for (g in 1:", n_groups, ") {"),
              paste0("    ", u, "[g] <- ", u_std, "[g] / sqrt(", tau_u, ")"),
              paste0("  }")
            )
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
          }

          # Handle multi-tree with generic structure naming
          s_name <- if (length(structure_names) > 0) {
            structure_names[1]
          } else {
            "struct"
          }

          s_lvl <- get_struct_lvl(s_name, hierarchical_info)
          s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
          s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)

          prec_name <- paste0("Prec_", s_name)
          prec_index <- if (is_multi_structure) {
            paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, ", K]")
          } else {
            paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, "]")
          }

          # Sum all terms: overdispersion + random effects + MAG terms
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) {
            paste0(" + ", paste(mag_terms, collapse = " + "))
          } else {
            ""
          }

          loop_bound <- get_loop_bound(response, hierarchical_info)

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for Poisson: ", response),
            paste0("  for (i in 1:", loop_bound, ") {"),
            paste0("    ", epsilon, "[i] ~ dnorm(0, ", tau_e, ")"),
            paste0("    ", err, "[i] <- ", epsilon, "[i]", total_u, total_mag),
            paste0("  }")
          )
        } else {
          # Non-optimised Poisson (Marginal or Latent)
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) {
            paste0(" + ", paste(mag_terms, collapse = " + "))
          } else {
            ""
          }

          if (independent) {
            # handled above, but just in case
          } else {
            # Standard GLMM definition of err
            model_lines <- c(
              model_lines,
              paste0("  # Non-optimised error for Poisson: ", response),
              paste0("  ", err, "[1:N] ~ dmnorm(zero_vec[], ", tau, ")"),
              paste0("  for (i in 1:N) {"),
              paste0(
                "    err_combined_",
                response,
                suffix,
                "[i] <- ",
                err,
                "[i]",
                total_mag
              ),
              paste0("  }")
            )
            # Need to update the likelihood logic to use err_combined if needed,
            # but for now let's just definition err[i] directly if possible.
            # Actually, err[i] ~ dmnorm is a vector.
          }
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

        # Initialize random effects accumulator for psi
        total_u <- ""
        # Initialize structure effect accumulator
        total_u <- ""
        if (optimise && !is.null(structures)) {
          for (s_name in names(structures)) {
            s_suffix <- paste0("_", s_name) # Always use structure name suffix
            u_var <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            # Structure Generic
            def <- jags_structure_definition(
              structures[[s_name]],
              variable_name = u_var,
              optimize = optimise,
              precision_parameter = tau_u
            )

            if (!is.null(def$error_prior)) {
              model_lines <- c(model_lines, def$error_prior)
            }
            total_u <- paste0(total_u, " + ", u_var, "[i]")
          }
        }

        # Start Observation Loop and Link Psi
        model_lines <- c(
          model_lines,
          paste0("  # Occupancy Model for ", response),
          paste0("  for (i in 1:N) {"),
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
          model_lines <- c(model_lines, def$model_code)
        } else {
          stop(paste("Unknown distribution or missing module for:", dist))
        }

        model_lines <- c(model_lines, "  }")
      } else if (dist == "negbinomial" || dist == "zinb") {
        # Negative Binomial error term: err[1:N]
        # Single phylogenetic effect (like Poisson/ordinal)
        err <- paste0("err_", response, suffix)

        if (optimise) {
          # Optimized Random Effects
          epsilon <- paste0("epsilon_", response, suffix)
          tau_e <- paste0("tau_e_", response, suffix)

          total_u <- ""

          # Phylogenetic / N-dim Structures
          for (s_name in structure_names) {
            if (
              !is_valid_structure_mapping(
                get_struct_lvl(s_name, hierarchical_info),
                get_var_level(response, hierarchical_info),
                hierarchical_info
              )
            ) {
              next
            }
            # Structure properties
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_suffix <- paste0("_", s_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            prec_name <- paste0("Prec_", s_name)
            prec_index <- if (multi.tree && is_multi_structure) {
              paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, ", K]")
            } else {
              paste0(prec_name, "[1:", s_bound, ", 1:", s_bound, "]")
            }

            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                u_std,
                "[1:",
                s_bound,
                "] ~ dmnorm(",
                s_zeros,
                "[1:",
                s_bound,
                "], ",
                prec_index,
                ")"
              ),
              paste0(
                "  for (i in 1:N) { ",
                u,
                "[i] <- ",
                u_std,
                "[",
                get_struct_index(s_name, response, hierarchical_info),
                "] / sqrt(",
                tau_u,
                ") }"
              )
            )
            total_u <- paste0(total_u, " + ", u, "[i]")
          }

          # Random Group Structures
          for (r_name in random_structure_names) {
            s_suffix <- paste0("_", r_name)
            u_std <- paste0("u_std_", response, suffix, s_suffix)
            u <- paste0("u_", response, suffix, s_suffix)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            n_groups <- paste0("N_", r_name)

            # Smart Group Index Logic for Hierarchy
            group_idx <- paste0("group_", r_name, "[i]")
            if (!is.null(hierarchical_info)) {
              resp_bound <- get_loop_bound(response, hierarchical_info)
              grp_bound <- get_loop_bound(r_name, hierarchical_info)

              # Strip N_ prefix
              resp_lvl <- sub("^N_", "", resp_bound)
              grp_lvl <- sub("^N_", "", grp_bound)

              # Resolve "N" to finest level name
              finest_lvl <- trimws(strsplit(hierarchical_info$hierarchy, ">")[[
                1
              ]])
              finest_lvl <- finest_lvl[length(finest_lvl)]

              if (resp_lvl == "N") {
                resp_lvl <- finest_lvl
              }
              if (grp_lvl == "N") {
                grp_lvl <- finest_lvl
              }

              # If levels differ, use nested transitive index
              # This creates u[ group_var[ bridge_idx[i] ] ]
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
            prec_name <- paste0("Prec_", r_name)
            zeros_name <- paste0("zeros_", r_name)

            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                u_std,
                "[1:",
                n_groups,
                "] ~ dmnorm(",
                zeros_name,
                "[1:",
                n_groups,
                "], ",
                prec_name,
                "[1:",
                n_groups,
                ", 1:",
                n_groups,
                "])"
              ),
              paste0("  for (g in 1:", n_groups, ") {"),
              paste0("    ", u, "[g] <- ", u_std, "[g] / sqrt(", tau_u, ")"),
              paste0("  }")
            )
            total_u <- paste0(total_u, " + ", u, "[", group_idx, "]")
          }

          # For Negative Binomial / ZINB, we do NOT add an independent residual error (epsilon)
          # because the distribution already handles overdispersion via 'r' (size).
          # Adding epsilon creates identifiability issues (double overdispersion).

          model_lines <- c(
            model_lines,
            paste0("  # Random effects for Negative Binomial: ", response),
            paste0("  for (i in 1:N) {")
          )

          # Sum all terms: random effects + MAG terms
          mag_terms <- vars_error_terms[[response]]
          total_mag <- if (length(mag_terms) > 0) {
            paste0(" + ", paste(mag_terms, collapse = " + "))
          } else {
            ""
          }

          if (total_u == "" && total_mag == "") {
            model_lines <- c(model_lines, paste0("    ", err, "[i] <- 0"))
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

    # Define Structure Error (if structure exists)
    phylo_term <- ""
    if (length(structure_names) > 0) {
      # Use first structure for induced correlations
      s_name <- structure_names[1]
      prec_name <- paste0("Prec_", s_name)
      err_phylo <- paste0("err_", s_name, "_", var)

      if (optimise) {
        # Optimised: use Prec_<structure> (for MAG/induced correlations)
        prec_index <- if (is_multi_structure) {
          paste0(prec_name, "[1:", loop_bound, ", 1:", loop_bound, ", K]")
        } else {
          paste0(prec_name, "[1:", loop_bound, ", 1:", loop_bound, "]")
        }

        model_lines <- c(
          model_lines,
          paste0(
            "  ",
            get_precision_prior(paste0("tau_", s_name, "_", var), var)
          ),
          paste0(
            "  ",
            err_phylo,
            "[1:",
            loop_bound,
            "] ~ dmnorm(zero_vec[1:", # Assume zero_vec is large enough or make new one?
            # Actually zero_vec is created as 1:N.
            # If loop_bound < N, zero_vec[1:loop_bound] is valid.
            loop_bound,
            "], ",
            "tau_",
            s_name,
            "_",
            var,
            " * ",
            prec_index,
            ")"
          )
        )
      } else {
        # Non-optimised: use Prec_<structure> directly
        prec_index <- if (is_multi_structure) {
          paste0(prec_name, "[1:", loop_bound, ", 1:", loop_bound, ", K]")
        } else {
          paste0(prec_name, "[1:", loop_bound, ", 1:", loop_bound, "]")
        }

        model_lines <- c(
          model_lines,
          paste0(
            "  ",
            get_precision_prior(paste0("tau_", s_name, "_", var), var)
          ),
          paste0(
            "  ",
            err_phylo,
            "[1:",
            loop_bound,
            "] ~ dmnorm(zero_vec[1:",
            loop_bound,
            "], ",
            "tau_",
            s_name,
            "_",
            var,
            " * ",
            prec_index,
            ")"
          )
        )
      }
      phylo_term <- paste0(" + ", err_phylo, "[i]")
    }

    # Define Observation Precision (Estimated)
    # This allows proper variance decomposition between:
    # - Structural effects (beta coefficients)
    # - Correlation structure (err_res from Wishart)
    # - Observation noise (tau_obs)
    model_lines <- c(
      model_lines,
      get_precision_prior(paste0("tau_obs_", var), var),
      paste0("  sigma_obs_", var, " <- 1/sqrt(tau_obs_", var, ")")
    )

    # Build sum string
    sum_res_errs <- paste(err_terms, collapse = " + ")

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
        phylo_term,
        " + ",
        sum_res_errs,
        ", tau_obs_",
        var,
        ")"
      )
    )
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
        phylo_term,
        " + ",
        sum_res_errs,
        ", tau_obs_",
        var,
        ")"
      )
    )
    model_lines <- c(model_lines, paste0("  }"))
  }

  # Measurement error / Variability
  if (length(variability_list) > 0) {
    model_lines <- c(model_lines, "  # Measurement error / Variability")

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
          paste0("  for (i in 1:N) {"),
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
        model_lines <- c(model_lines, paste0("  }"))
      } else if (type == "reps") {
        # For repeated measures, we need log_lik for EACH observation
        # Then sum them per individual for WAIC
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {")
        )
        # Replaced imperative accumulation with declarative summation
        model_lines <- c(
          model_lines,
          paste0("    for (j in 1:N_reps_", var, "[i]) {"),
          paste0(
            "      ",
            var,
            "_obs[i, j] ~ dnorm(",
            var,
            "[i], ",
            var,
            "_tau)"
          ),
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
          ),
          paste0("    }"),
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
        model_lines <- c(
          model_lines,
          paste0("  }"),
          get_precision_prior(paste0(var, "_tau"), var),
          paste0("  ", var, "_sigma <- 1/sqrt(", var, "_tau)")
        )
      } else {
        warning(paste("Unknown variability type:", type, "for variable", var))
      }
    }
  }

  model_lines <- c(model_lines, "  # Priors for structural parameters")

  # Priors for alpha, lambda, tau, sigma
  for (response in names(response_counter)) {
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
      default_alpha <- "dnorm(0, 1.0E-6)"
      if (dist %in% logit_dists || is_occupancy_aux) {
        default_alpha <- "dnorm(0, 1)"
      }
      model_lines <- c(
        model_lines,
        paste0("  ", get_prior(alpha_name, default_alpha))
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
        optimise &&
        (length(structure_names) > 0 || length(random_structure_names) > 0) &&
        dist != "occupancy"

      if (
        (is.null(vars_with_na) || !response %in% vars_with_na || optimise) &&
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
                "  tau_e_",
                response,
                suffix,
                " <- ",
                prec,
                " # Fixed residual variance"
              )
            )
          } else {
            model_lines <- c(
              model_lines,
              paste0(
                "  ",
                get_precision_prior(
                  paste0("tau_e_", response, suffix),
                  response
                )
              )
            )
          }
          model_lines <- c(
            model_lines,
            paste0(
              "  sigma",
              response,
              suffix,
              " <- 1/sqrt(tau_e_",
              response,
              suffix,
              ")"
            )
          )
        } else if (optimise) {
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
                "  tau_e_",
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
              model_lines <- c(
                model_lines,
                paste0(
                  "  ",
                  get_precision_prior(
                    paste0("tau_e_", response, suffix),
                    response
                  )
                )
              )
            }
          }

          # Generate Structure Priors
          for (s_name in structure_names) {
            if (
              !is_valid_structure_mapping(
                get_struct_lvl(s_name, hierarchical_info),
                get_var_level(response, hierarchical_info),
                hierarchical_info
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

            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            model_lines <- c(
              model_lines,
              paste0("  ", get_precision_prior(tau_u, response)),
              paste0(
                "  sigma_",
                response,
                suffix,
                s_suffix,
                " <- 1/sqrt(",
                tau_u,
                ")"
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
            model_lines <- c(
              model_lines,
              paste0(
                "  lambda",
                response,
                suffix,
                " <- (1/",
                tau_u_name,
                ") / ((1/",
                tau_u_name,
                ") + (1/tau_e_",
                response,
                suffix,
                "))"
              )
            )
          }

          # Random Effects Priors
          for (r_name in random_structure_names) {
            s_suffix <- paste0("_", r_name)
            tau_u <- paste0("tau_u_", response, suffix, s_suffix)

            model_lines <- c(
              model_lines,
              paste0("  ", get_precision_prior(tau_u, response)),
              paste0(
                "  sigma_",
                response,
                suffix,
                s_suffix,
                " <- 1/sqrt(",
                tau_u,
                ")"
              )
            )
          }

          # Residual Sigma (only if tau_e exists)
          if (!dist %in% c("negbinomial", "zinb")) {
            model_lines <- c(
              model_lines,
              paste0(
                "  sigma_",
                response,
                "_res <- 1/sqrt(tau_e_",
                response,
                suffix,
                ")"
              )
            )
          }
        } else {
          # Marginal Priors (lambda, tau)
          model_lines <- c(
            model_lines,
            paste0(
              "  ",
              get_prior(paste0("lambda", response, suffix), "dunif(0, 1)")
            ),
            paste0(
              "  ",
              get_precision_prior(paste0("tau", response, suffix), response)
            ),
            paste0(
              "  sigma",
              response,
              suffix,
              " <- 1/sqrt(tau",
              response,
              suffix,
              ")"
            )
          )
        }
      }
    }
  }

  # Priors for multinomial parameters (arrays)
  for (response in names(response_counter)) {
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
            "    tau_e_",
            response,
            "[k] <- ",
            prec,
            " # Fixed"
          )
        } else {
          tau_line <- paste0(
            "    ",
            paste0(
              "    ",
              get_precision_prior(paste0("tau_e_", response, "[k]"), response)
            )
          )
        }

        model_lines <- c(
          model_lines,
          paste0("  # Independent Priors for ", response, " (Multinomial)"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 1.0E-6)"),
          tau_line,
          "  }"
        )
      } else if (optimise) {
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
            "    tau_e_",
            response,
            "[k] <- ",
            prec,
            " # Fixed"
          )
        } else {
          tau_line <- paste0(
            "    ",
            get_precision_prior(paste0("tau_e_", response, "[k]"), response)
          )
        }

        model_lines <- c(
          model_lines,
          # Priors for ", response, " (Multinomial)"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 1.0E-6)"),

          # Residual error
          tau_line,

          # Random effects priors (loop over structures)
          {
            prior_lines <- c()

            # Phylogenetic / N-dim Structures
            for (s_name in structure_names) {
              if (
                !is_valid_structure_mapping(
                  get_struct_lvl(s_name, hierarchical_info),
                  get_var_level(response, hierarchical_info),
                  hierarchical_info
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
              tau_u <- paste0("tau_u_", response, s_suffix) # Suffix handled outside loop
              # Wait, suffix is from k=1..response_counter. But here we are in k=2..K_var (categories).
              # The external loop is 'response', but 'response' is unique per equation group.
              # Multinomial doesn't support repeats in this logic currently (response_counter[[response]] is likely 1).
              # We use [k] for category index.

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
            for (r_name in random_structure_names) {
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
              "[k]) + (1/tau_e_",
              response,
              "[k]))"
            )
          } else {
            NULL
          },
          "  }"
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0("  # Priors for ", response, " (Multinomial)"),
          paste0("  for (k in 2:", K_var, ") {"),
          paste0("    alpha_", response, "[k] ~ dnorm(0, 1.0E-6)"),
          paste0("    lambda_", response, "[k] ~ dunif(0, 1)"),
          paste0(
            "    ",
            paste0(
              "    ",
              get_precision_prior(paste0("tau_", response, "[k]"), response)
            )
          ),
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
              paste0("  for (k in 2:", K_var, ") {"),
              paste0(
                "    ",
                get_prior(paste0(beta_name, "[k]"), "dnorm(0, 1.0E-6)")
              ),
              "  }"
            )
          }
        }
      }
    }
  }

  # Priors for ordinal parameters (cutpoints + variance components)
  for (response in names(response_counter)) {
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
                  paste0("tau_e_", response, suffix),
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
                  paste0("tau_e_", response, suffix),
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
                ") + (1/tau_e_",
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
              paste0("  ", get_prior(beta_name, "dnorm(0, 1.0E-6)"))
            )
          }
        }
      }
    }
  }

  # Priors for Negative Binomial parameters (size only - others handled in main loop)
  for (response in names(response_counter)) {
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
      paste0("  alpha_", var, " ~ dnorm(0, 1.0E-6)")
    )
  }

  # Priors for Zero-Inflation parameters (psi)
  for (response in names(response_counter)) {
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

  for (response in names(response_counter)) {
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

  for (beta in unique_betas) {
    # Isolate response part of beta name: beta_Response_Predictor
    # Assuming names are strictly beta_Resp_Pred, but Resp could contain underscores if user did so.
    # A safer way is using param_map, but we are in a simple loop here.
    # Let's try to match against names(dist_list).

    default_beta <- "dnorm(0, 1.0E-6)"

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

    model_lines <- c(
      model_lines,
      paste0("  ", get_prior(beta, default_beta))
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
        model_lines <- c(model_lines, paste0("  ", line))
      }
    }
  }

  # Phylogenetic tree selection (if multi.tree)
  if (multi.tree) {
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
      # Only use legacy GLMM blocking if optimisation implies marginal approach (i.e. optimise=FALSE)
      use_glmm <- (!is.null(vars_with_na) &&
        response %in% vars_with_na &&
        !optimise)

      if (dist == "gaussian" && !use_glmm) {
        if (!optimise) {
          if (multi.tree) {
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
        tau_struct <- paste0("tau_", s_name, "_", response, suffix)
        tau_res <- paste0("tau_res_", response, suffix)

        # Priors
        model_lines <- c(
          model_lines,
          paste0("  ", get_precision_prior(tau_res, var)),
          paste0("  ", get_precision_prior(tau_struct, var)),
          # Calculate lambda for reporting
          paste0(
            "  lambda",
            response,
            suffix,
            " <- (1/",
            tau_struct,
            ") / ((1/",
            tau_struct,
            ") + (1/",
            tau_res,
            "))"
          )
        )

        if (is_multi_structure) {
          model_lines <- c(
            model_lines,
            paste0(
              "  TAU_",
              s_name,
              "_",
              response,
              suffix,
              " <- ",
              tau_struct,
              " * inverse(multiVCV[,,K])"
            ),
            paste0(
              "  ",
              err,
              "[1:N] ~ dmnorm(",
              mu_err,
              "[], TAU_",
              s_name,
              "_",
              response,
              suffix,
              ")"
            )
          )
        } else {
          model_lines <- c(
            model_lines,
            paste0(
              "  TAU_",
              s_name,
              "_",
              response,
              suffix,
              " <- ",
              tau_struct,
              " * inverse(VCV)"
            ),
            paste0(
              "  ",
              err,
              "[1:N] ~ dmnorm(",
              mu_err,
              "[], TAU_",
              s_name,
              "_",
              response,
              suffix,
              ")"
            )
          )
        }
      } else if (dist == "binomial") {
        # GLMM covariance for binomial error term
        if (!optimise) {
          # Marginal Formulation
          if (multi.tree) {
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
        # tau_u priors for optimise=TRUE are already generated in the main priors block (lines ~2293-2311)
        # No need to duplicate here.
      } else if (dist == "occupancy") {
        # Occupancy priors are handled by the generic loop above (lines 2243+ for structures)
        # No specific residual variance (tau_e) needed.
      } else if (dist == "multinomial") {
        # Multinomial covariance
        # We need TAU[,,k] for each k
        K_var <- paste0("K_", response)

        if (!optimise) {
          # k=1 is reference category (fixed to identity)
          # k>=2 have estimated phylogenetic signal
          if (multi.tree) {
            model_lines <- c(
              model_lines,
              paste0("  # Covariance matrices for multinomial"),
              "  # Reference category k=1",
              "  for (i in 1:N) {",
              "    for (j in 1:N) {",
              paste0("      Mlam_", response, "[i,j,1] <- ID[i,j]"),
              "    }",
              "  }",
              paste0(
                "  TAU_",
                tolower(response),
                "_",
                suffix,
                "[1:N,1:N,1] <- ID[1:N,1:N]"
              ),
              "  # Estimated categories k>=2",
              paste0("  for (k in 2:", K_var, ") {"),
              "    for (i in 1:N) {",
              "      for (j in 1:N) {",
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
                "[1:N,1:N,k] <- inverse(Mlam_",
                response,
                "[,,k])"
              ),
              "  }"
            )
          } else {
            model_lines <- c(
              model_lines,
              paste0("  # Covariance matrices for multinomial"),
              "  # Reference category k=1",
              paste0(
                "  TAU_",
                tolower(response),
                "_",
                suffix,
                "[1:N,1:N,1] <- ID[1:N,1:N]"
              ),
              "  # Estimated categories k>=2",
              paste0("  for (k in 2:", K_var, ") {"),
              paste0(
                "    TAU_",
                tolower(response),
                "_",
                suffix,
                "[1:N, 1:N, k] <- inverse(lambda_",
                response,
                "[k] * VCV + (1 - lambda_",
                response,
                "[k]) * ID)"
              ),
              "  }"
            )
          }
        }
      }
    }
  }

  # Covariance for correlated vars (phylogenetic part)
  # Only use VCV approach when optimise=FALSE; optimised models use eigendecomposition
  if (!is.null(induced_correlations) && !optimise) {
    s_name <- if (length(structure_names) > 0) structure_names[1] else "struct"
    for (var in correlated_vars) {
      if (is_multi_structure) {
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
    model_lines <- c(model_lines, "  # Predictor priors for imputation")
  }
  non_response_vars <- setdiff(all_vars, names(response_counter))
  for (var in non_response_vars) {
    # Check if this is a latent variable
    is_latent <- !is.null(latent) && var %in% latent

    # Skip fully observed predictors (not latent and no missing data)
    # Only generate imputation for: latent variables OR variables with missing data
    if (!is_latent && (is.null(vars_with_na) || !var %in% vars_with_na)) {
      next # Skip this predictor - it's fully observed data
    }

    model_lines <- c(
      model_lines,
      paste0("  for (i in 1:", get_loop_bound(var, hierarchical_info), ") {"),
      paste0("    mu", var, "[i] <- 0"),
      paste0("  }")
    )

    if (independent) {
      # Independent imputation (i.i.d normal)
      if (is_latent && standardize_latent) {
        # Use N(0,1) prior for standardized latent variable
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
            "[i] ~ dnorm(0, 1)  # Standardized latent variable"
          ),
          paste0("  }")
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0(
            "  for (i in 1:",
            get_loop_bound(var, hierarchical_info),
            ") {"
          ),
          paste0("    ", var, "[i] ~ dnorm(mu", var, "[i], tau_e_", var, ")"),
          paste0("  }")
        )
      }
    } else if (optimise) {
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
            "[i] ~ dnorm(0, 1)  # Standardized latent variable"
          ),
          paste0("  }")
        )
      } else {
        # Standard random effects formulation
        additive_terms <- ""

        if (!is.null(structures)) {
          for (s_name in names(structures)) {
            s_suffix <- paste0("_", s_name)
            u_var <- paste0("u_", var, s_suffix)
            tau_u <- paste0("tau_u_", var, s_suffix)

            # Call Generic to define the error term u ~ dmnorm(...)
            def <- jags_structure_definition(
              structures[[s_name]],
              variable_name = u_var,
              optimize = optimise,
              precision_parameter = tau_u
            )

            if (!is.null(def$error_prior)) {
              model_lines <- c(model_lines, def$error_prior)
            }

            additive_terms <- paste0(additive_terms, " + ", u_var, "[i]")
          }
        }

        tau_e <- paste0("tau_e_", var)

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
    } else {
      # Standard MVN (Marginal)
      # If latent variable with standardize_latent = TRUE, use simple N(0,1) prior
      if (is_latent && standardize_latent) {
        model_lines <- c(
          model_lines,
          paste0("  for (i in 1:N) {"),
          paste0(
            "    ",
            var,
            "[i] ~ dnorm(0, 1)  # Standardized latent variable"
          ),
          paste0("  }")
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0(
            "  ",
            var,
            "[1:N] ~ dmnorm(mu",
            var,
            "[1:N], TAU",
            tolower(var),
            ")"
          )
        )
      }
    }

    # For latent variables, fix tau = 1 (standardize)
    # For observed predictors, estimate tau
    # Skip if standardize_latent = TRUE (variance already set in N(0,1) prior)
    if (is_latent && !standardize_latent) {
      if (optimise) {
        model_lines <- c(
          model_lines,
          paste0("  # Latent variable: standardized (var = 1)"),
          paste0("  lambda", var, " ~ dunif(0, 1)"),
          # In random effects: Var = (1/tau_u)*V + (1/tau_e)*I
          # We want Var = lambda*V + (1-lambda)*I
          # So 1/tau_u = lambda => tau_u = 1/lambda
          # And 1/tau_e = 1-lambda => tau_e = 1/(1-lambda)
          paste0("  tau_u_", var, " <- 1/lambda", var),
          paste0("  tau_e_", var, " <- 1/(1-lambda", var, ")"),
          paste0("  sigma", var, " <- 1") # Fixed to 1
        )
      } else {
        model_lines <- c(
          model_lines,
          paste0("  # Latent variable: standardized (var = 1)"),
          paste0("  lambda", var, " ~ dunif(0, 1)"),
          paste0("  tau", var, " <- 1  # Fixed for identification"),
          paste0("  sigma", var, " <- 1/sqrt(tau", var, ")")
        )
      }
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
            "  tau_e_",
            var,
            " <- ",
            prec,
            " # Fixed residual variance"
          )
        } else {
          tau_line <- paste0(
            "  ",
            get_precision_prior(paste0("tau_e_", var), var)
          )
        }

        # Independent Predictor Prior
        model_lines <- c(
          model_lines,
          tau_line,
          paste0("  sigma", var, " <- 1/sqrt(tau_e_", var, ")")
        )
      } else if (optimise) {
        if (TRUE) {
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
              "  tau_e_",
              var,
              " <- ",
              prec,
              " # Fixed residual variance"
            )
          } else {
            tau_line <- paste0(
              "  ",
              get_precision_prior(paste0("tau_e_", var), var)
            )
          }

          # Multiple Structures: Estimate independent variance components
          if (needs_residual_variance) {
            model_lines <- c(
              model_lines,
              tau_line
            )
          }

          for (s_name in structure_names) {
            if (
              !is_valid_structure_mapping(
                get_struct_lvl(s_name, hierarchical_info),
                get_var_level(response, hierarchical_info),
                hierarchical_info
              )
            ) {
              next
            }
            # Structure properties
            s_lvl <- get_struct_lvl(s_name, hierarchical_info)
            s_bound <- if (is.null(s_lvl)) "N" else paste0("N_", s_lvl)
            s_zeros <- if (is.null(s_lvl)) "zeros" else paste0("zeros_", s_lvl)
            s_suffix <- paste0("_", s_name)
            tau_u <- paste0("tau_u_", var, s_suffix)
            model_lines <- c(
              model_lines,
              paste0("  ", get_precision_prior(tau_u, var)),
              paste0("  sigma_", var, s_suffix, " <- 1/sqrt(", tau_u, ")")
            )
          }

          if (needs_residual_variance) {
            model_lines <- c(
              model_lines,
              paste0("  sigma_", var, "_res <- 1/sqrt(tau_e_", var, ")")
            )
          }
        } else {
          # Single Structure (Legacy behavior with lambda partitioning)
          if (!needs_residual_variance) {
            # No residual variance (e.g. occupancy), so just estimate tau_u directly
            model_lines <- c(
              model_lines,
              paste0("  ", get_precision_prior(paste0("tau_u_", var), var)),
              paste0("  sigma", var, " <- 1/sqrt(tau_u_", var, ")")
            )
          } else {
            # Normal Gaussian case: Partition total variance tau with lambda
            model_lines <- c(
              model_lines,
              paste0("  lambda", var, " ~ dunif(0, 1)"),
              paste0("  ", get_precision_prior(paste0("tau", var), var)),
              paste0("  tau_u_", var, " <- tau", var, "/lambda", var),
              paste0("  tau_e_", var, " <- tau", var, "/(1-lambda", var, ")"),
              paste0("  sigma", var, " <- 1/sqrt(tau", var, ")")
            )
          }
        }
      } else {
        # Single Structure Case for Imputed Variable
        s_name <- structure_names[1]
        s_suffix <- paste0("_", s_name)
        tau_u <- paste0("tau_u_", var, s_suffix)

        # Need to define tau_u for this variable
        model_lines <- c(
          model_lines,
          paste0("  ", get_precision_prior(tau_u, var)),
          paste0("  sigma_", var, s_suffix, " <- 1/sqrt(", tau_u, ")")
        )

        # Also generate tau_e if needed (same logic as above)
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
          prec <- 1 / val
          model_lines <- c(
            model_lines,
            paste0("  tau_e_", var, " <- ", prec, " # Fixed")
          )
        } else {
          model_lines <- c(
            model_lines,
            paste0("  ", get_precision_prior(paste0("tau_e_", var), var))
          )
        }

        # Calculate lambda for compatibility
        model_lines <- c(
          model_lines,
          paste0(
            "  lambda",
            var,
            s_suffix,
            " <- (1/",
            tau_u,
            ") / ((1/",
            tau_u,
            ") + (1/tau_e_",
            var,
            "))"
          )
        )
      }
    }

    # Only generate TAU matrix with VCV when there is a structure defined
    if (!optimise && length(structure_names) > 0) {
      if (multi.tree) {
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

  model_lines <- c(model_lines, "}")
  model_string <- paste(model_lines, collapse = "\n")

  # Convert param_map to data frame
  param_map_df <- do.call(
    rbind,
    lapply(param_map, as.data.frame, stringsAsFactors = FALSE)
  )

  return(list(model = model_string, parameter_map = param_map_df))
}
