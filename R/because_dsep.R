#' Extract d-separation statements from a structural equation model
#'
#' This function takes a set of structural equations defining a causal model
#' and returns the conditional independence statements (d-separation or m-separation tests)
#' implied by the model structure. If latent variables are specified, the function
#' uses the MAG (Mixed Acyclic Graph) approach by Shipley and Douma (2021)
#' to account for unmeasured latent variables.
#'
#' @param equations A list of model formulas (one per structural equation),
#'   e.g., \code{list(Y ~ X1 + X2, Z ~ Y)}.
#' @param latent Optional character vector of latent (unmeasured) variable names.
#'   If provided, the function converts the DAG to a MAG and returns m-separation tests.
#'
#' @return If \code{latent} is NULL, returns a list of formulas representing
#'   conditional independence tests. If \code{latent} is specified, returns a list with:
#'   \itemize{
#'     \item \code{tests}: List of m-separation test formulas
#'     \item \code{correlations}: List of variable pairs with induced correlations
#'   }
#'
#' @details
#' The function implements the basis set approach to d-separation testing
#' (Shipley 2000, 2009, 2016). For standard DAGs without latent variables, it identifies
#' pairs of non-adjacent variables and creates conditional independence tests.
#'
#' When latent variables are specified, the function uses the DAG-to-MAG conversion
#' (Shipley & Douma 2021) to identify m-separation statements and induced correlations
#' among observed variables that arise from shared latent common causes.
#'
#' Deterministic nodes (interaction terms such as \code{A:B}, and arithmetic
#' transformations such as \code{I(A^2)}) are kept as **explicit intermediate
#' nodes** in the DAG, following the D-separation extension of
#' Geiger, Verma & Pearl (1990).  This ensures that the basis set includes
#' independence tests that condition on the deterministic term itself
#' (e.g. \eqn{TL \perp BM \mid \{BM{:}M\}}), which would be silently dropped
#' if the interaction were collapsed to its component variables.
#'
#' @references
#' Geiger, D., Verma, T., & Pearl, J. (1990). Identifying independence in
#' Bayesian Networks. \emph{Networks}, 20(5), 507--534.
#'
#' Shipley, B. (2000). A new inferential test for path models based on
#' directed acyclic graphs. Structural Equation Modeling, 7(2), 206-218.
#'
#' Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
#' context. Ecology, 90(2), 363-368.
#'
#' Shipley, B. (2016). Cause and Correlation in Biology (2nd ed.).
#' Cambridge University Press.
#'
#' Shipley, B., & Douma, J. C. (2021). Testing Piecewise Structural Equations
#' Models in the Presence of Latent Variables and Including Correlated Errors.
#' Structural Equation Modeling: A Multidisciplinary Journal, 28(4), 582-589.
#' https://doi.org/10.1080/10705511.2020.1871355
#'

#'
#' @examples
#' # Standard DAG
#' equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)
#' ind_tests <- because_dsep(equations)
#'
#' # With latent variable
#' equations_latent <- list(X ~ Quality, Y ~ Quality)
#' result <- because_dsep(equations_latent, latent = "Quality")
#' # result$tests: m-separation tests
#' # result$correlations: induced correlation between X and Y
#' @param poly_terms Internal list of polynomial terms.
#' @param categorical_vars Character vector of categorical variable names.
#' @param family Named character vector of family/distribution for response variables.
#' @param quiet Logical; if FALSE (default), print the basis set and MAG structure.
#'   If TRUE, suppress informational output.
#' @param random_terms Optional list of random effects (group, type) parsed from equations.
#' @param hierarchical_info Internal argument used to pass data hierarchy information
#'   (levels, grouping variables) for multiscale d-separation
#'   tests (following Shipley 2009).
#' @export
#' @importFrom stats formula terms as.formula
because_dsep <- function(
  equations,
  latent = NULL,
  random_terms = list(),
  hierarchical_info = NULL,
  poly_terms = NULL,
  categorical_vars = NULL,
  family = NULL,
  quiet = FALSE
) {
  # Support passing a fitted 'because' object directly
  if (inherits(equations, "because")) {
    obj <- equations
    equations <- obj$equations
    if (is.null(latent) && !is.null(obj$latent)) latent <- obj$latent
    if (is.null(hierarchical_info) && !is.null(obj$hierarchical_info)) hierarchical_info <- obj$hierarchical_info
  }

  # If no latents, use standard DAG d-separation
  if (is.null(latent)) {
    return(dsep_standard(
      equations,
      random_terms = random_terms,
      hierarchical_info = hierarchical_info,
      poly_terms = poly_terms,
      categorical_vars = categorical_vars,
      family = family,
      quiet = quiet
    ))
  }

  # With latents: use MAG m-separation
  return(dsep_with_latents(
    equations,
    latent,
    random_terms = random_terms,
    hierarchical_info = hierarchical_info,
    poly_terms = poly_terms,
    categorical_vars = categorical_vars,
    family = family,
    quiet = quiet
  ))
}


# --- Hierarchical Helpers for D-Sep ---

# Get the level name for a variable
get_var_level_dsep <- function(var, hierarchical_info) {
  if (is.null(hierarchical_info) || is.null(hierarchical_info$levels)) {
    return(NULL)
  }
  for (lvl in names(hierarchical_info$levels)) {
    if (var %in% hierarchical_info$levels[[lvl]]) {
      return(lvl)
    }
  }
  return(NULL)
}

# Determine if two hierarchy levels belong to separate, non-nested branches.
# The hierarchy string has the form "a > b > c; d > e" meaning the semicolon
# separates completely independent chains. Two levels are "orthogonal" if they
# appear in different chains and neither is an ancestor of the other.
are_levels_orthogonal <- function(lvl_a, lvl_b, hierarchy_str) {
  if (is.null(lvl_a) || is.null(lvl_b) || lvl_a == lvl_b) return(FALSE)
  paths <- strsplit(hierarchy_str, "\\s*;\\s*")[[1]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1]]))
  # Find which paths each level belongs to
  path_a <- which(sapply(paths, function(p) lvl_a %in% p))
  path_b <- which(sapply(paths, function(p) lvl_b %in% p))
  # Orthogonal = they are in different chains
  if (length(path_a) == 0 || length(path_b) == 0) return(FALSE)
  !any(path_a %in% path_b)
}

#' Check if a d-sep test is a cross-hierarchy test (trivially satisfied)
#'
#' A test Response _||_ FocalPredictor | \code{\{...\}} is "cross-hierarchy" when the
#' response and the focal predictor live in orthogonal hierarchical branches
#' (e.g., one is species-level, the other is site-level) AND the conditioning
#' set does not contain a variable from a level that connects the two branches
#' (typically the observation level). Such tests are trivially independent by
#' design of the hierarchical data structure and cannot be run in JAGS because
#' no cross-level index exists between orthogonal branches.
#'
#' @param test_eq A formula with attribute "test_var" naming the focal predictor.
#' @param hierarchical_info List with $levels (named list) and $hierarchy (string).
#' @return TRUE if the test should be skipped; FALSE otherwise.
#' @keywords internal
is_cross_hierarchy_test <- function(test_eq, hierarchical_info) {
  if (is.null(hierarchical_info) ||
      is.null(hierarchical_info$levels) ||
      is.null(hierarchical_info$hierarchy)) {
    return(FALSE)
  }
  test_var <- attr(test_eq, "test_var")
  
  # FALLBACK: If attribute is missing, assume first predictor on RHS is the focal one
  if (is.null(test_var)) {
    rhs_vars <- all.vars(test_eq)[-1]
    if (length(rhs_vars) > 0) test_var <- rhs_vars[1]
  }

  if (is.null(test_var)) return(FALSE)

  # Response is the LHS
  resp <- as.character(test_eq)[2]
  resp <- sub("^psi_", "", resp)

  resp_lvl <- get_var_level_dsep(resp, hierarchical_info)
  pred_lvl <- get_var_level_dsep(test_var, hierarchical_info)

  if (is.null(resp_lvl) || is.null(pred_lvl)) return(FALSE)
  if (resp_lvl == pred_lvl) return(FALSE)

  # Are the two levels in separate hierarchy chains?
  if (!are_levels_orthogonal(resp_lvl, pred_lvl, hierarchical_info$hierarchy)) {
    return(FALSE)
  }

  # Even if orthogonal, check whether any conditioning variable bridges them
  # via the obs (or another shared) level. If the conditioning set contains a
  # variable at the obs level (or any level that appears in BOTH chains), the
  # path could be opened. For safety, only skip if conditioning set has NO
  # bridging variable.
  rhs <- as.character(test_eq)[3]
  cond_terms <- trimws(strsplit(rhs, "\\+")[[1]])
  # Strip random effects (1 | X)
  cond_terms <- cond_terms[!grepl("^\\s*1\\s*\\|", cond_terms)]
  # Strip I(...) interaction terms (they are derived, not level-assigning)
  cond_terms <- cond_terms[!grepl("^\\s*I\\(", cond_terms)]
  cond_terms <- trimws(gsub("\\(.*\\)", "", cond_terms))
  cond_terms <- cond_terms[nchar(cond_terms) > 0]

  # Check if any conditioning variable is at a level that appears in both chains
  hierarchy_str <- hierarchical_info$hierarchy
  paths <- strsplit(hierarchy_str, "\\s*;\\s*")[[1]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1]]))
  path_a <- which(sapply(paths, function(p) resp_lvl %in% p))
  path_b <- which(sapply(paths, function(p) pred_lvl %in% p))

  for (cv in cond_terms) {
    cv_lvl <- get_var_level_dsep(cv, hierarchical_info)
    if (is.null(cv_lvl)) next
    cv_in_a <- any(sapply(paths[path_a], function(p) cv_lvl %in% p))
    cv_in_b <- any(sapply(paths[path_b], function(p) cv_lvl %in% p))

    # A variable bridges only if it appears in BOTH lineages. 
    # Example: 'obs' appears in both 'site > survey > obs' and 'species > obs'.
    # A variable at 'site' level does not bridge 'survey' to 'species'.
    if (cv_in_a && cv_in_b) return(FALSE)
  }

  return(TRUE)
}

#' Get inherited random terms for a variable based on hierarchy
#'
#' @param var Variable name
#' @param hierarchical_info List with $levels, $hierarchy, and $link_vars
#' @return List of random terms (group, type)
#' @keywords internal
get_inherited_random_terms <- function(var, hierarchical_info) {
  if (is.null(hierarchical_info) ||
      is.null(hierarchical_info$levels) ||
      is.null(hierarchical_info$hierarchy) ||
      is.null(hierarchical_info$link_vars)) {
    return(list())
  }

  lvl <- get_var_level_dsep(var, hierarchical_info)
  if (is.null(lvl)) return(list())

  # Find all ancestors of this level in the hierarchy string
  # (Splitting by ; handles separate chains)
  paths <- strsplit(hierarchical_info$hierarchy, "\\s*;\\s*")[[1]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1]]))

  ancestors <- character(0)
  for (path in paths) {
    idx <- match(lvl, path)
    if (!is.na(idx) && idx > 1) {
      ancestors <- c(ancestors, path[1:(idx - 1)])
    }
  }
  ancestors <- unique(ancestors)

  inherited <- list()
  for (anc in ancestors) {
    if (anc %in% names(hierarchical_info$link_vars)) {
      group_var <- hierarchical_info$link_vars[[anc]]
      if (!is.null(group_var)) {
        inherited[[length(inherited) + 1]] <- list(
          group = group_var,
          type = "intercept"
        )
      }
    }
  }

  return(inherited)
}

# Check if a structure's level can map to a response's level
is_valid_structure_mapping_dsep <- function(s_lvl, r_lvl, hierarchical_info) {
  if (is.null(s_lvl) || is.null(r_lvl) || s_lvl == r_lvl) {
    return(TRUE)
  }
  if (is.null(hierarchical_info) || is.null(hierarchical_info$hierarchy)) {
    return(TRUE)
  }

  paths <- strsplit(hierarchical_info$hierarchy, "\\s*;\\s*")[[1]]
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

# Standard d-separation for DAGs (using dagitty)
dsep_standard <- function(
  equations,
  random_terms = list(),
  hierarchical_info = NULL,
  poly_terms = NULL,
  categorical_vars = NULL,
  family = NULL,
  quiet = FALSE
) {
  # Extract grouping variables from random terms to exclude from DAG
  grouping_vars <- NULL
  if (length(random_terms) > 0) {
    grouping_vars <- unique(sapply(random_terms, function(x) x$group))
  }

  # Extract ALL deterministic terms (interactions + I() calls).
  # These are kept as explicit intermediate nodes in the DAG, following
  # Geiger, Verma & Pearl (1990), so that tests like TL ⊥ BM | {BM:M}
  # appear in the basis set.
  if (!is.null(poly_terms)) {
    all_poly_terms <- poly_terms
  } else {
    all_poly_terms <- get_all_polynomial_terms(equations)
  }
  all_det_terms <- extract_deterministic_terms(equations)

  # Grouping variables are excluded from the DAG; deterministic nodes are NOT
  # excluded — they are kept as explicit nodes (see equations_to_dag).
  exclude_vars <- grouping_vars

  # Normalize equations for DAG using Modular Family Logic
  norm_equations <- equations
  if (!is.null(family)) {
    unique_dists <- unique(family)
    for (dist in unique_dists) {
      fam_obj <- get_family_object(dist)
      vars_with_dist <- names(family)[family == dist]
      for (v in vars_with_dist) {
        norm_equations <- transform_graph_for_dsep(
          fam_obj,
          norm_equations,
          variable = v
        )
      }
    }
  }

  dag <- equations_to_dag(
    norm_equations,
    exclude_vars = exclude_vars,
    deterministic_terms = all_det_terms
  )

  # Use dagitty to get the correct d-separation basis set
  dag_str <- dag_matrix_to_dagitty(dag)
  d_obj <- dagitty::dagitty(dag_str)
  # Compute true Shipley (2000) minimal basis set
  sorted <- dagitty::topologicalOrdering(d_obj)
  sorted_nodes <- names(sorted)[order(unlist(sorted))]
  
  basis <- list()
  for (i in seq_along(sorted_nodes)) {
    v_i <- sorted_nodes[i]
    parents_vi <- dagitty::parents(d_obj, v_i)
    
    if (i > 1) {
      predecessors <- sorted_nodes[1:(i-1)]
      non_parents <- setdiff(predecessors, parents_vi)
      
      for (v_j in non_parents) {
        parents_vj <- dagitty::parents(d_obj, v_j)
        basis[[length(basis) + 1]] <- c(v_i, v_j, unique(c(parents_vi, parents_vj)))
      }
    }
  }
  # Build combined exclusion list:
  #   1. Deterministic interaction/I() nodes from extract_deterministic_terms.
  #   2. Poly term internal names (e.g. age_pow2) — needed when because() has
  #      already expanded equations (I(age^2) -> age_pow2) before calling here,
  #      so extract_deterministic_terms no longer finds them.
  excl_names <- unique(c(
    if (!is.null(all_det_terms) && length(all_det_terms) > 0) {
      sapply(all_det_terms, function(x) x$internal_name)
    } else {
      character(0)
    },
    if (!is.null(all_poly_terms) && length(all_poly_terms) > 0) {
      sapply(all_poly_terms, function(x) x$internal_name)
    } else {
      character(0)
    }
  ))
  if (length(excl_names) > 0 && !is.null(basis)) {
    basis <- Filter(
      function(test) {
        !(test[1] %in% excl_names || test[2] %in% excl_names)
      },
      basis
    )
  }

  # Polynomial-term injection into conditioning sets (for I(x^2) type terms)
  # With explicit det nodes in the DAG this is now largely handled structurally,
  # but we keep the injection for any residual polynomial terms.
  if (!is.null(all_poly_terms) && length(basis) > 0) {
    basis <- lapply(basis, function(test) {
      if (length(test) > 2) {
        cond_vars <- test[3:length(test)]
        new_cond_vars <- cond_vars
        for (cv in cond_vars) {
          for (pt in all_poly_terms) {
            if (pt$base_var == cv) {
              new_cond_vars <- c(new_cond_vars, pt$internal_name)
            }
          }
        }
        return(c(test[1:2], unique(new_cond_vars)))
      } else {
        return(test)
      }
    })
  }

  # Filter out random effect grouping variables from conditioning sets
  if (length(random_terms) > 0 && !is.null(basis)) {
    basis <- lapply(basis, function(test) {
      if (length(test) > 2) {
        cond_vars <- test[3:length(test)]
        filtered_cond <- cond_vars[!cond_vars %in% grouping_vars]
        return(c(test[1:2], filtered_cond))
      }
      return(test)
    })
  }

  # Convert basis set to formula list
  # Root nodes (no parents in the DAG) should always be predictors, not responses.
  root_vars <- rownames(dag)[colSums(dag) == 0]

  tests <- mag_basis_to_formulas(
    basis,
    categorical_vars = categorical_vars,
    family = family,
    deterministic_terms = all_det_terms,
    root_vars = root_vars,
    hierarchical_info = hierarchical_info,
    d_obj = d_obj,
    quiet = quiet
  )

  # Append random terms if relevant
  if (length(random_terms) > 0 && length(tests) > 0) {
    new_tests <- list()
    for (t_idx in seq_along(tests)) {
      t_eq <- tests[[t_idx]]
      resp <- as.character(t_eq)[2]

      base_resp <- sub("^psi_", "", resp)
      vocab_rand <- Filter(
        function(x) {
          # Must match response name
          if (x$response != resp && x$response != base_resp) {
            return(FALSE)
          }

          # If hierarchical info is present, check level compatibility
          if (!is.null(hierarchical_info)) {
            r_lvl <- get_var_level_dsep(base_resp, hierarchical_info)
            g_lvl <- get_var_level_dsep(x$group, hierarchical_info)

            # If both levels are known, check compatibility (group must be coarser or equal)
            if (!is.null(r_lvl) && !is.null(g_lvl)) {
              if (
                !is_valid_structure_mapping_dsep(
                  g_lvl,
                  r_lvl,
                  hierarchical_info
                )
              ) {
                return(FALSE)
              }
            }
          }
          return(TRUE)
        },
        random_terms
      )

      # --- Automated Hierarchical Random Effects ---
      # If hierarchical info is present, automatically add random effects for
      # grouping variables at ancestor levels to avoid pseudo-replication.
      # e.g., if response is at 'survey' level, add '(1 | Site)'.
      inherited_rand <- get_inherited_random_terms(base_resp, hierarchical_info)

      # Combine and deduplicate by group name
      combined_rand <- vocab_rand
      for (ir in inherited_rand) {
        if (!any(sapply(combined_rand, function(x) x$group == ir$group))) {
          # Only add if the grouping variable is NOT the predictor itself
          test_var <- attr(t_eq, "test_var")
          if (is.null(test_var) || ir$group != test_var) {
            combined_rand[[length(combined_rand) + 1]] <- ir
          }
        }
      }
      vocab_rand <- combined_rand

      if (length(vocab_rand) > 0) {
        rand_str <- paste(
          sapply(vocab_rand, function(rt) {
            paste0("(1 | ", rt$group, ")")
          }),
          collapse = " + "
        )

        f_str <- paste(deparse(t_eq), collapse = " ")
        f_str <- paste0(f_str, " + ", rand_str)
        new_eq <- as.formula(f_str)
        attr(new_eq, "test_var") <- attr(t_eq, "test_var")
        new_tests[[t_idx]] <- new_eq
      } else {
        new_tests[[t_idx]] <- t_eq
      }
    }
    tests <- new_tests
  }

  # Print basis set if not quiet
  if (!quiet) {
    cat("Basis Set for DAG:", "\n")
    cat(
      "I(X,Y|Z) means X is d-separated from Y given the set Z in the DAG",
      "\n"
    )
    if (length(tests) == 0) {
      cat("No elements in the basis set", "\n")
    } else {
      for (test in tests) {
        # [FIX] Pass random_terms to show grouping variables in conditioning set
        cat(format_dsep_test(test, random_terms = random_terms), "\n")
      }
    }
  }

  return(tests)
}

# M-separation for MAGs (with latent variables)
dsep_with_latents <- function(
  equations,
  latent,
  random_terms = list(),
  hierarchical_info = NULL,
  poly_terms = NULL,
  categorical_vars = NULL,
  family = NULL,
  quiet = FALSE
) {
  # --- Modular Graph Transformation ---
  augmented_equations <- equations
  if (!is.null(family)) {
    unique_dists <- unique(family)
    for (dist in unique_dists) {
      fam_obj <- get_family_object(dist)
      vars_with_dist <- names(family)[family == dist]
      for (v in vars_with_dist) {
        augmented_equations <- transform_graph_for_dsep(
          fam_obj,
          augmented_equations,
          variable = v
        )
      }
    }
  }

  # Extract grouping variables from random terms to exclude from DAG
  grouping_vars <- NULL
  if (length(random_terms) > 0) {
    grouping_vars <- unique(sapply(random_terms, function(x) x$group))
  }

  # Extract ALL deterministic terms (interactions + I() calls).
  # Kept as explicit nodes in the DAG following Geiger, Verma & Pearl (1990).
  if (!is.null(poly_terms)) {
    all_poly_terms <- poly_terms
  } else {
    all_poly_terms <- get_all_polynomial_terms(equations)
  }
  all_det_terms <- extract_deterministic_terms(equations)

  # Only grouping vars are excluded; deterministic nodes remain in the graph.
  exclude_vars <- grouping_vars

  dag <- equations_to_dag(
    augmented_equations,
    exclude_vars = exclude_vars,
    deterministic_terms = all_det_terms
  )

  # Use dagitty to handle latent variables natively
  dag_str <- dag_matrix_to_dagitty(dag)
  d_obj <- dagitty::dagitty(dag_str)
  if (!is.null(latent)) {
    valid_latents <- intersect(latent, names(d_obj))
    if (length(valid_latents) > 0) {
      dagitty::latents(d_obj) <- valid_latents
    }
  }
  
  # --- Compute optimal MAG basis set (Shipley & Douma 2021) ---
  # Find which observed nodes share a latent parent
  latent_sharers <- list()
  if (!is.null(latent) && length(latent) > 0) {
    for (l in latent) {
      if (l %in% names(d_obj)) {
        children <- dagitty::children(d_obj, l)
        for (c in children) {
          if (is.null(latent_sharers[[c]])) {
            latent_sharers[[c]] <- character(0)
          }
          latent_sharers[[c]] <- unique(c(latent_sharers[[c]], setdiff(children, c)))
        }
      }
    }
  }
  
  # Extract all possible implied independencies (dagitty natively marginalizes latents)
  indeps <- unclass(dagitty::impliedConditionalIndependencies(d_obj))
  
  # Build combined exclusion list for deterministic/poly terms
  excl_names <- unique(c(
    if (!is.null(all_det_terms) && length(all_det_terms) > 0) {
      sapply(all_det_terms, function(x) x$internal_name)
    } else {
      character(0)
    },
    if (!is.null(all_poly_terms) && length(all_poly_terms) > 0) {
      sapply(all_poly_terms, function(x) x$internal_name)
    } else {
      character(0)
    }
  ))

  # Group by pair and apply collinearity penalty
  pairs <- list()
  for (i in seq_along(indeps)) {
    x <- as.character(unlist(indeps[[i]]$X))
    y <- as.character(unlist(indeps[[i]]$Y))
    z <- as.character(unlist(indeps[[i]]$Z))
    
    # Filter out deterministic/polynomial nodes from focal pairs
    if (x %in% excl_names || y %in% excl_names) next
    
    xy <- sort(c(x, y))
    key <- paste(xy, collapse = "_||_")
    
    if (is.null(pairs[[key]])) {
      pairs[[key]] <- list()
    }
    
    col_pen <- 0
    if (length(z) > 0) {
      sharers_x <- if (!is.null(latent_sharers[[x]])) latent_sharers[[x]] else character(0)
      sharers_y <- if (!is.null(latent_sharers[[y]])) latent_sharers[[y]] else character(0)
      all_sharers <- unique(c(sharers_x, sharers_y))
      col_pen <- sum(z %in% all_sharers)
    }
    size_pen <- length(z)
    
    pairs[[key]][[length(pairs[[key]]) + 1]] <- list(
      X = x, Y = y, Z = z, col_pen = col_pen, size_pen = size_pen
    )
  }
  
  basis <- list()
  for (key in names(pairs)) {
    options <- pairs[[key]]
    options <- options[order(sapply(options, function(x) x$col_pen), sapply(options, function(x) x$size_pen))]
    best <- options[[1]]
    basis[[length(basis) + 1]] <- c(best$X, best$Y, best$Z)
  }

  # Polynomial-term injection in conditioning sets
  if (!is.null(all_poly_terms) && length(basis) > 0) {
    basis <- lapply(basis, function(test) {
      if (length(test) > 2) {
        cond_vars <- test[3:length(test)]
        new_cond_vars <- cond_vars
        for (cv in cond_vars) {
          for (pt in all_poly_terms) {
            if (pt$base_var == cv) {
              new_cond_vars <- c(new_cond_vars, pt$internal_name)
            }
          }
        }
        return(c(test[1:2], unique(new_cond_vars)))
      } else {
        return(test)
      }
    })
  }

  # Filter out random effect grouping variables from conditioning sets
  if (length(random_terms) > 0 && !is.null(basis)) {
    grouping_vars <- unique(sapply(random_terms, function(x) x$group))
    basis <- lapply(basis, function(test) {
      if (length(test) > 2) {
        cond_vars <- test[3:length(test)]
        filtered_cond <- cond_vars[!cond_vars %in% grouping_vars]
        c(test[1:2], filtered_cond)
      } else {
        test
      }
    })
  }

  # Identify variables that are direct children of latent variables
  latent_children <- character(0)
  if (!is.null(latent) && length(latent) > 0) {
    all_vars <- rownames(dag)
    for (lat in latent) {
      if (lat %in% all_vars) {
        children <- all_vars[dag[lat, ] == 1]
        latent_children <- unique(c(latent_children, children))
      }
    }
  }

  # Convert to formula format
  # Root nodes (no parents in the DAG) should always be predictors, not responses.
  root_vars <- rownames(dag)[colSums(dag) == 0]

  tests <- mag_basis_to_formulas(
    basis,
    latent_children = latent_children,
    categorical_vars = categorical_vars,
    family = family,
    deterministic_terms = all_det_terms,
    root_vars = root_vars,
    hierarchical_info = hierarchical_info,
    latent_sharers = latent_sharers,
    d_obj = d_obj,
    quiet = quiet
  )

  # Save tests without random effects for clean display

  # Save tests without random effects for clean display
  tests_for_display <- tests

  # Append random terms to MAG tests
  if (length(random_terms) > 0 && length(tests) > 0) {
    new_tests <- list()
    for (t_idx in seq_along(tests)) {
      t_eq <- tests[[t_idx]]
      resp <- as.character(t_eq)[2]

      base_resp <- sub("^psi_", "", resp)
      vocab_rand <- Filter(
        function(x) {
          # Must match response name
          if (x$response != resp && x$response != base_resp) {
            return(FALSE)
          }

          # If hierarchical info is present, check level compatibility
          if (!is.null(hierarchical_info)) {
            r_lvl <- get_var_level_dsep(base_resp, hierarchical_info)
            g_lvl <- get_var_level_dsep(x$group, hierarchical_info)

            if (!is.null(r_lvl) && !is.null(g_lvl)) {
              if (
                !is_valid_structure_mapping_dsep(
                  g_lvl,
                  r_lvl,
                  hierarchical_info
                )
              ) {
                return(FALSE)
              }
            }
          }
          return(TRUE)
        },
        random_terms
      )

      # --- Automated Hierarchical Random Effects ---
      inherited_rand <- get_inherited_random_terms(base_resp, hierarchical_info)

      # Combine and deduplicate by group name
      combined_rand <- vocab_rand
      for (ir in inherited_rand) {
        if (!any(sapply(combined_rand, function(x) x$group == ir$group))) {
          # Only add if the grouping variable is NOT the predictor itself
          test_var <- attr(t_eq, "test_var")
          if (is.null(test_var) || ir$group != test_var) {
            combined_rand[[length(combined_rand) + 1]] <- ir
          }
        }
      }
      vocab_rand <- combined_rand

      if (length(vocab_rand) > 0) {
        rand_str <- paste(
          sapply(vocab_rand, function(rt) {
            paste0("(1 | ", rt$group, ")")
          }),
          collapse = " + "
        )
        f_str <- paste(deparse(t_eq), collapse = " ")
        f_str <- paste0(f_str, " + ", rand_str)
        new_eq <- as.formula(f_str)
        new_tests[[t_idx]] <- new_eq
      } else {
        new_tests[[t_idx]] <- t_eq
      }
    }
    tests <- new_tests
  }

  # Extract bidirected edges (induced correlations)
  # Induced correlations in dagitty are handled via bidirected edges in the MAG
  mag_obj <- dagitty::toMAG(d_obj)
  mag_edges <- dagitty::edges(mag_obj)
  correlations <- list()
  if (nrow(mag_edges) > 0) {
    bidirected <- mag_edges[mag_edges$e == "<->", ]
    if (nrow(bidirected) > 0) {
      for (i in seq_len(nrow(bidirected))) {
        correlations[[length(correlations) + 1]] <- c(bidirected$v[i], bidirected$w[i])
      }
    }
  }

  # Print basis set if not quiet
  if (!quiet) {
    cat("Basis Set for MAG:", "\n")
    cat(
      "I(X,Y|Z) means X is m-separated from Y given the set Z in the MAG",
      "\n"
    )
    if (length(tests_for_display) == 0) {
      cat("No elements in the basis set", "\n")
    } else {
      for (test in tests_for_display) {
        # [FIX] Pass random_terms to show grouping variables in conditioning set
        cat(format_dsep_test(test, random_terms = random_terms), "\n")
      }
    }
  }

  return(list(
    tests = tests,
    correlations = correlations,
    mag = d_obj
  ))
}

# Helper to format a d-sep test for printing
format_dsep_test <- function(test, random_terms = NULL) {
  # Extract response and test variable
  resp <- as.character(test)[2]
  test_var <- attr(test, "test_var")
  
  # Use deparse to get full RHS
  rhs_str <- paste(deparse(test[[3]]), collapse = " ")
  # Split by + and trim
  rhs_parts <- trimws(strsplit(rhs_str, "\\+")[[1]])
  
  # The conditioning set is everything in the RHS except the test_var
  cond_set <- rhs_parts[rhs_parts != test_var]
  
  if (is.null(test_var)) {
     # Fallback for old tests without attribute
     vars <- all.vars(test)
     test_var <- setdiff(vars, resp)[1]
     cond_set <- setdiff(vars, c(resp, test_var))
  }
  
  if (length(cond_set) == 0) {
    return(paste0("I( ", resp, " , ", test_var, " |  )"))
  } else {
    # Sort for canonical representation (fixed terms first, then random)
    is_rand <- grepl("\\|", cond_set)
    fixed_cond <- sort(cond_set[!is_rand])
    rand_cond <- sort(cond_set[is_rand])
    
    full_cond <- paste(c(fixed_cond, rand_cond), collapse = ", ")
    return(paste0("I( ", resp, " , ", test_var, " | ", full_cond, " )"))
  }
}


#' plot_dsep
#'
#' Creates a caterpillar plot (point and whisker) of the regression coefficients
#' from all d-separation tests. A horizontal red line at zero helps visually
#' assess which independence claims are fulfilled (95% CI includes zero) or
#' violated (95% CI excludes zero).
#'
#' @param object A `because` object fitted with \code{dsep = TRUE}.
#' @param ... Additional arguments.
#' @param prob Numeric; probability mass for the credibility interval (default 0.95).
#'
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' # Plot results for a fitted model
#' plot_dsep(fit)
#' }
#' @rdname plot_dsep
#' @export
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_hline coord_flip labs theme_minimal theme
plot_dsep.because <- function(object, prob = 0.95, ...) {
  has_dsep <- !is.null(object$dsep) &&
    (isTRUE(object$dsep) ||
       (is.list(object$dsep) && !is.null(object$dsep$results)))
  if (!has_dsep) {
    stop("plot_dsep requires a 'because' object fitted with dsep = TRUE.")
  }

  # summary.because handles the complex renaming and dummy variable matching
  s <- summary(object, prob = prob)

  if (is.null(s$results) || nrow(s$results) == 0) {
    stop("No d-separation test results found in model object.")
  }

  res <- s$results

  # Ensure ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_dsep. Please install it.")
  }

  # For multinomial responses, multiple parameters share the same Test string
  # (one per category). Disambiguate by appending the specific Parameter name
  # so that each category gets its own row in the caterpillar plot.
  dup_tests <- duplicated(res$Test) | duplicated(res$Test, fromLast = TRUE)
  if (any(dup_tests)) {
    res$Test[dup_tests] <- paste0(res$Test[dup_tests], " (", res$Parameter[dup_tests], ")")
  }

  # Create Simplified labels for the plot: Response _||_ TestVar
  # This strips the conditioning set (e.g. | {X,Y,Z}) for better readability.
  # We match the specific formatting used in summary.because.
  res$Label <- gsub(" \\| \\{.*?\\}", "", res$Test)
  res$Label <- trimws(res$Label)

  # Symmetric axis: always show both sides of zero for fair visual comparison.
  # Use coord_flip(ylim=) rather than scale_y_continuous(limits=) to zoom
  # the coordinate system without clipping whiskers at the scale level.
  max_abs <- max(abs(c(res$LowerCI, res$UpperCI, res$Estimate)), na.rm = TRUE)
  axis_lim <- c(-max_abs, max_abs) * 1.05  # 5% padding

  # Create Plot
  p <- ggplot2::ggplot(
    res,
    ggplot2::aes(
      x = stats::reorder(Label, seq_len(nrow(res))),
      y = Estimate,
      ymin = LowerCI,
      ymax = UpperCI
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey50",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::geom_pointrange(linewidth = 0.8, size = 0.5, fatten = 3) +
    ggplot2::coord_flip(ylim = axis_lim) +
    ggplot2::labs(
      title = "d-separation Independence Tests",
      subtitle = paste0("Caterpillar plot of path coefficients with ", prob * 100, "% Bayesian Credibility Intervals"),
      x = "Conditional Independence Claim",
      y = "Estimated Beta (Effect Size)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 11),
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}


# ==========================================================================
# Cross-Scale D-Sep: PGLS routing
# ==========================================================================
# When the focal predictor lives at a *coarser* hierarchical level than the
# response (e.g., species-level Body_Mass_s tested against obs-level Abundance),
# a standard observation-level GLMM is non-identifiable: the group-level fixed
# effect and the group-level random effect compete for the same variance.
#
# To avoid MCMC convergence issues, we run the test at the predictor's own scale:
#   1. Aggregate the response to the predictor's level.
#   2. Include only conditioning variables at or above that level.
#   3. Fit GLS with the appropriate correlation structure (phylo / spatial).
# ==========================================================================

#' Detect whether a d-sep test crosses hierarchical scales in the same chain
#'
#' A test is "cross-scale" when the focal predictor is at a coarser level than
#' the response within the SAME hierarchical chain (not orthogonal branches).
#' E.g., species-level Body_Mass_s predicting obs-level Abundance in the chain
#' "site > survey > obs; species > obs".
#'
#' @param test_eq A formula carrying the attribute \code{test_var}.
#' @param hierarchical_info List with \code{$levels}, \code{$hierarchy},
#'   \code{$link_vars}.
#' @return A list. If cross-scale: \code{$is_crossscale = TRUE} plus
#'   \code{$response}, \code{$response_level}, \code{$predictor_level},
#'   \code{$test_var}. Otherwise \code{$is_crossscale = FALSE}.
#' @keywords internal
detect_crossscale_dsep <- function(test_eq, hierarchical_info) {
  if (is.null(hierarchical_info) ||
      is.null(hierarchical_info$levels) ||
      is.null(hierarchical_info$hierarchy)) {
    return(list(is_crossscale = FALSE))
  }

  test_var <- attr(test_eq, "test_var")
  if (is.null(test_var)) {
    rhs_vars <- all.vars(test_eq)[-1]
    if (length(rhs_vars) > 0L) test_var <- rhs_vars[1L]
  }
  if (is.null(test_var)) return(list(is_crossscale = FALSE))

  resp     <- as.character(test_eq)[2L]
  resp     <- sub("^psi_", "", resp)

  resp_lvl <- get_var_level_dsep(resp,     hierarchical_info)
  pred_lvl <- get_var_level_dsep(test_var, hierarchical_info)

  if (is.null(resp_lvl) || is.null(pred_lvl)) return(list(is_crossscale = FALSE))
  if (resp_lvl == pred_lvl)                   return(list(is_crossscale = FALSE))

  # Parse hierarchy chains (semicolon-separated, ">" = coarser > finer)
  paths <- strsplit(hierarchical_info$hierarchy, "\\s*;\\s*")[[1L]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1L]]))

  for (path in paths) {
    resp_idx <- match(resp_lvl, path)
    pred_idx <- match(pred_lvl, path)
    if (!is.na(resp_idx) && !is.na(pred_idx) && pred_idx < resp_idx) {
      # Predictor is coarser (smaller index) than response in the same chain
      return(list(
        is_crossscale   = TRUE,
        response        = resp,
        response_level  = resp_lvl,
        predictor_level = pred_lvl,
        test_var        = test_var
      ))
    }
  }
  list(is_crossscale = FALSE)
}

.cs_flatten_all <- function(original_data) {
  if (is.data.frame(original_data)) return(original_data)
  if (!is.list(original_data)) return(original_data)
  
  row_counts <- sapply(original_data, function(x) if (is.data.frame(x)) nrow(x) else 0)
  base_idx <- which.max(row_counts)
  if (length(base_idx) == 0) return(NULL)
  
  base_tbl <- original_data[[base_idx]]
  merged_tables <- c(base_idx)
  
  repeat {
    merged_in_this_pass <- FALSE
    for (i in seq_along(original_data)) {
      if (i %in% merged_tables) next
      tbl <- original_data[[i]]
      if (!is.data.frame(tbl)) next
      
      shared <- intersect(names(base_tbl), names(tbl))
      if (length(shared) > 0) {
        new_cols <- setdiff(names(tbl), names(base_tbl))
        if (length(new_cols) > 0) {
          base_tbl <- merge(base_tbl, tbl[, c(shared, new_cols), drop = FALSE], by = shared, all.x = TRUE)
        }
        merged_tables <- c(merged_tables, i)
        merged_in_this_pass <- TRUE
      }
    }
    # Count dataframes
    n_dfs <- sum(sapply(original_data, is.data.frame))
    if (!merged_in_this_pass || length(merged_tables) == n_dfs) break
  }
  
  base_tbl
}

#' Run a cross-scale d-sep test via aggregation and PGLS/GLS
#'
#' Implements the scale-aware d-sep testing strategy for hierarchical models:
#' \enumerate{
#'   \item Aggregate the response variable to the predictor's hierarchical level
#'         (log-mean for Poisson; arithmetic mean for Gaussian).
#'   \item Select conditioning variables at or above the predictor's level.
#'   \item Fit a GLS with \code{ape::corPagel} (Pagel's lambda estimated by ML)
#'         if a phylogenetic tree matching the predictor's level is found in
#'         \code{structure}; otherwise fall back to ordinary OLS.
#'   \item Return synthetic MCMC samples (normal approximation to the GLS
#'         posterior of the focal coefficient) so that the result slots directly
#'         into because()'s existing downstream summary / plot code.
#' }
#'
#' @param i            Integer index of the d-sep test.
#' @param test_eq      Formula carrying \code{attr(., "test_var")}.
#' @param cs_info      List returned by \code{detect_crossscale_dsep()}.
#' @param original_data Hierarchical data list or flat data.frame.
#' @param hierarchical_info List with \code{$levels}, \code{$hierarchy},
#'   \code{$link_vars}.
#' @param structure    Named list of structure objects (phylo trees, spatial
#'   matrices), as passed to \code{because()}.
#' @param family       Named character vector mapping response names to families.
#' @param n.iter       Synthetic draws per MCMC chain (default 1000).
#' @param n.chains     Number of synthetic chains (default 3).
#' @param quiet        Suppress informational messages.
#' @return List with \code{$samples} (coda::mcmc.list), \code{$param_map},
#'   \code{$model}, \code{$test_index}.
#' @keywords internal
run_crossscale_dsep_pgls <- function(
  i,
  test_eq,
  cs_info,
  original_data,
  hierarchical_info,
  structure = NULL,
  family    = NULL,
  engine    = "numpyro",
  n.iter    = 1000L,
  n.chains  = 3L,
  quiet     = FALSE
) {
  resp       <- cs_info$response
  test_var   <- cs_info$test_var
  pred_level <- cs_info$predictor_level
  link_vars  <- hierarchical_info$link_vars
  group_col  <- link_vars[[pred_level]]

  if (is.null(group_col)) {
    stop(sprintf(
      "Cross-scale d-sep: no link_var defined for level '%s'. Cannot aggregate.",
      pred_level
    ))
  }

  # ---- 1. Use the globally flattened dataset ----
  flat <- .cs_flatten_all(original_data)

  # ---- 2. Aggregate response to the predictor's group level ----
  y_raw   <- flat[[resp]]
  grp     <- as.character(flat[[group_col]])
  fam_str <- if (!is.null(family) && resp %in% names(family)) family[[resp]] else "gaussian"

  agg_y <- if (fam_str %in% c("poisson", "negbin", "zip", "zinb")) {
    tapply(log(pmax(y_raw, 0) + 0.5), grp, mean, na.rm = TRUE)
  } else if (fam_str %in% c("lognormal", "gamma", "exponential")) {
    tapply(log(pmax(y_raw, 1e-6)), grp, mean, na.rm = TRUE)
  } else if (fam_str %in% c("binomial", "bernoulli", "beta")) {
    tapply(y_raw, grp, function(x) {
      p <- mean(x, na.rm = TRUE)
      p <- pmax(0.01, pmin(0.99, p)) # Bound to prevent log(0)
      log(p / (1 - p)) # Logit transform
    })
  } else {
    tapply(y_raw, grp, mean, na.rm = TRUE)
  }
  group_ids <- names(agg_y)

  agg_df <- data.frame(
    stringsAsFactors = FALSE,
    setNames(list(group_ids, as.numeric(agg_y)), c(group_col, ".response"))
  )

  # ---- 3. Aggregate all conditioning variables to the group level ----
  # Exclude random effect grouping variables (which are just the link_vars)
  rhs_vars <- setdiff(all.vars(test_eq)[-1L], c(group_col, unlist(link_vars)))
  
  for (v in rhs_vars) {
    if (v %in% names(flat)) {
      v_raw <- flat[[v]]
      agg_v <- tapply(v_raw, grp, function(x) {
        if (is.numeric(x)) {
          mean(x, na.rm = TRUE)
        } else {
          # For categorical predictors, take the most frequent value (mode)
          tbl <- sort(table(x), decreasing = TRUE)
          if (length(tbl) > 0) names(tbl)[1] else NA
        }
      })
      # Ensure order matches agg_df
      agg_df[[v]] <- as.vector(agg_v[agg_df[[group_col]]])
    }
  }

  agg_df <- agg_df[stats::complete.cases(agg_df), , drop = FALSE]

  # ---- 4. Reformat aggregated dataset for flat Bayesian run ----
  avail_preds <- intersect(rhs_vars, names(agg_df))
  if (length(avail_preds) == 0L) {
    stop("No predictors available in aggregated data for cross-scale PGLS.")
  }
  
  # Ensure the response has its original name for the because() call
  names(agg_df)[names(agg_df) == ".response"] <- resp
  
  # For the phylogeny mapping in because(), if the data is flat, the row names or 
  # an identifier column matching the tree tips must be present.
  # The identifier column is group_col, so we rename it to the expected raw format.
  # because() expects '.raw_Species' for species mapping. We use group_col dynamically.
  raw_group_col <- paste0(".raw_", group_col)
  names(agg_df)[names(agg_df) == group_col] <- raw_group_col
  
  # Ensure we only keep the necessary columns to avoid unrelated data issues
  keep_cols <- c(raw_group_col, resp, avail_preds)
  agg_df <- agg_df[, keep_cols, drop = FALSE]
  
  method_used <- sprintf("Bayesian PGLS/GLS (%s)", if (tolower(engine) == "jags") "JAGS" else if (tolower(engine) == "nimble") "NIMBLE" else "NumPyro")
  if (!quiet) {
    message(sprintf(
      "  [Cross-scale d-sep] '%s' aggregated to %s level (n=%d) -> %s ... ",
      resp, pred_level, nrow(agg_df), method_used
    ), appendLF = FALSE)
    flush.console()
  }

  # Build the flat equation
  gls_formula <- stats::as.formula(
    paste(resp, "~", paste(avail_preds, collapse = " + "))
  )

  # ---- 5. Prune structure to only what applies to this level ----
  pruned_structure <- list()
  if (!is.null(structure) && is.list(structure)) {
    groups_present <- as.character(agg_df[[raw_group_col]])
    for (s_name in names(structure)) {
      s_obj <- structure[[s_name]]
      # For phylo trees
      if (inherits(s_obj, "phylo")) {
        common <- intersect(s_obj$tip.label, groups_present)
        if (length(common) > 0) {
          pruned_structure[[s_name]] <- s_obj
        }
      } 
      # For spatial distance matrices
      else if (is.matrix(s_obj)) {
        common <- intersect(rownames(s_obj), groups_present)
        if (length(common) > 0) {
          pruned_structure[[s_name]] <- s_obj
        }
      }
    }
  }
  if (length(pruned_structure) == 0) pruned_structure <- NULL

  # ---- 6. Run native fully-Bayesian model via because() ----
  fit_bayesian <- tryCatch({
    because(
      equations = list(gls_formula),
      data = agg_df,
      structure = pruned_structure,  # Pass only applicable structures
      family = {
        if (!is.null(family)) {
          agg_fam <- as.list(family)
          agg_fam[[resp]] <- "gaussian"
          agg_fam
        } else {
          NULL
        }
      },
      dsep = FALSE,           # Just fit this single equation, don't generate more tests
      engine = engine,        # Use user-specified engine
      n.chains = n.chains,
      n.iter = n.iter,
      quiet = TRUE            # Suppress inner compilation messages
    )
  }, error = function(e) {
    warning(sprintf(
      "Cross-scale Bayesian PGLS failed (%s); returning empty samples for test %d.", conditionMessage(e), i
    ))
    NULL
  })

  param_name <- paste0("beta_", resp, "_", test_var)
  
  if (!quiet) {
    if (!is.null(fit_bayesian)) {
      message("Done.")
    } else {
      message("Failed.")
    }
  }

  # ---- 6. Extract MCMC chains and map parameters ----
  
  if (!is.null(fit_bayesian) && !is.null(fit_bayesian$samples)) {
    # Extract the chains for the focal predictor from the Bayesian run
    # because() names coefficients as `beta_Response_Predictor`
    if (inherits(fit_bayesian$samples, "mcmc.list")) {
      extracted_samples <- coda::mcmc.list(lapply(fit_bayesian$samples, function(x) {
        coda::mcmc(as.matrix(x)[, param_name, drop = FALSE])
      }))
    } else {
      extracted_samples <- fit_bayesian$samples[, param_name, drop = FALSE]
    }
  } else {
    # Fallback if the Bayesian model failed to compile/run
    extracted_samples <- lapply(seq_len(n.chains), function(ch) {
      mat <- matrix(NA_real_, nrow = n.iter, ncol = 1L, dimnames = list(NULL, param_name))
      coda::mcmc(mat)
    })
    extracted_samples <- coda::mcmc.list(extracted_samples)
  }

  # ---- 7. param_map (minimal columns consumed by downstream code) ----
  param_map <- data.frame(
    response       = resp,
    predictor      = test_var,
    parameter      = param_name,
    equation_index = i,
    type           = "beta",
    stringsAsFactors = FALSE
  )

  model_str <- sprintf(
    paste0(
      "# Cross-scale Bayesian PGLS d-sep test (test %d)\n",
      "# Response '%s' aggregated to '%s' level (n = %d)\n",
      "# Formula: %s ~ %s\n",
      "# Method: %s"
    ),
    i, resp, pred_level, nrow(agg_df),
    resp, paste(avail_preds, collapse = " + "),
    method_used
  )

  list(
    samples       = extracted_samples,
    param_map     = param_map,
    model         = model_str,
    test_index    = i,
    exact_p_value = NA_real_ # P-values are not standard in pure Bayesian outputs; credible intervals will be used
  )
}


# ---- Internal helpers (not exported) ------------------------------------

# Flatten a hierarchical data list to a single data.frame that contains both
# the response variable and the group column needed for aggregation.
.cs_flatten <- function(original_data, resp, group_col) {
  if (is.data.frame(original_data)) return(original_data)
  if (!is.list(original_data))
    stop("original_data must be a data.frame or a named list of data.frames.")

  # Find the table that holds the response
  base_tbl <- NULL
  for (tbl in original_data) {
    if (is.data.frame(tbl) && resp %in% names(tbl)) { base_tbl <- tbl; break }
  }
  if (is.null(base_tbl))
    stop(sprintf("Response variable '%s' not found in any data table.", resp))

  # If the group column is already present we are done
  if (group_col %in% names(base_tbl)) return(base_tbl)

  # Otherwise try to merge in the table that has group_col
  for (tbl in original_data) {
    if (!is.data.frame(tbl) || identical(tbl, base_tbl)) next
    if (group_col %in% names(tbl)) {
      shared   <- intersect(names(base_tbl), names(tbl))
      new_cols <- setdiff(names(tbl), names(base_tbl))
      if (length(shared) > 0L && length(new_cols) > 0L) {
        base_tbl <- merge(
          base_tbl,
          tbl[, c(shared, new_cols), drop = FALSE],
          by = shared, all.x = TRUE
        )
      }
      if (group_col %in% names(base_tbl)) break
    }
  }

  if (!group_col %in% names(base_tbl))
    stop(sprintf("Group column '%s' not found after attempting to flatten data.", group_col))

  base_tbl
}

# Find the data.frame (level table) that contains both group_col and test_var.
.cs_find_table <- function(original_data, group_col, test_var) {
  if (is.data.frame(original_data)) return(original_data)
  for (tbl in original_data) {
    if (is.data.frame(tbl) &&
        group_col %in% names(tbl) &&
        test_var  %in% names(tbl)) {
      return(tbl)
    }
  }
  stop(sprintf(
    "Cannot find a table containing both '%s' and '%s' for cross-scale PGLS.",
    group_col, test_var
  ))
}

# Build a phylogenetic correlation structure for nlme::gls.
# Returns list($cor, $agg_df) where $agg_df is row-aligned to tree tip order.
.cs_build_cor <- function(structure, agg_df, group_col) {
  if (is.null(structure) || !requireNamespace("ape", quietly = TRUE)) {
    return(list(cor = NULL, agg_df = agg_df))
  }

  groups <- as.character(agg_df[[group_col]])

  for (sn in names(structure)) {
    s_obj <- structure[[sn]]
    if (!inherits(s_obj, "phylo")) next

    tree   <- s_obj
    common <- intersect(tree$tip.label, groups)

    if (length(common) < 3L) {
      return(list(cor = NULL, agg_df = agg_df))
    }

    # Prune tree to the groups present in the data
    tree_p <- ape::drop.tip(tree, setdiff(tree$tip.label, groups))

    # Make ultrametric (required by corPagel)
    if (!ape::is.ultrametric(tree_p, tol = 1e-6)) {
      tree_p <- ape::compute.brlen(tree_p, method = "Grafen")
    }

    # Align agg_df rows to tree tip order
    rownames(agg_df) <- agg_df[[group_col]]
    agg_df_aligned   <- agg_df[tree_p$tip.label, , drop = FALSE]

    cor_struct <- ape::corPagel(1, phy = tree_p, form = stats::as.formula(paste("~", group_col)), fixed = FALSE)
    return(list(cor = cor_struct, agg_df = agg_df_aligned))
  }

  # No matching phylo structure found
  list(cor = NULL, agg_df = agg_df)
}
