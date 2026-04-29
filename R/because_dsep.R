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
    group_var <- hierarchical_info$link_vars[[anc]]
    if (!is.null(group_var)) {
      inherited[[length(inherited) + 1]] <- list(
        group = group_var,
        type = "intercept"
      )
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
  ici <- dagitty::impliedConditionalIndependencies(d_obj)
  
  # Convert dagitty output to the format expected by mag_basis_to_formulas
  basis <- lapply(ici, function(x) {
    c(as.character(x$X), as.character(x$Y), as.character(x$Z))
  })

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
    hierarchical_info = hierarchical_info
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
        cat(format_dsep_test(test), "\n")
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
    dagitty::latents(d_obj) <- latent
  }
  
  ici <- dagitty::impliedConditionalIndependencies(d_obj)
  basis <- lapply(ici, function(x) {
    c(as.character(x$X), as.character(x$Y), as.character(x$Z))
  })

  # Build combined exclusion list (same logic as dsep_standard)
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

  # --- Shipley & Douma (2021): Remove untestable basis set elements ---
  # Per the MAG approach, m-separation tests where a LATENT variable appears
  # in the CONDITIONING SET are not directly testable from observed data.
  # (You cannot condition on an unobserved variable in a regression.)
  # The latent's effect is already captured via bidirected edges in the MAG.
  # Filter these tests out so they are not run as JAGS sub-models.
  if (!is.null(latent) && length(latent) > 0 && !is.null(basis)) {
    n_before <- length(basis)
    basis <- Filter(
      function(test) {
        # Conditioning set is test[3:length(test)] (if present)
        if (length(test) > 2) {
          cond_vars <- test[3:length(test)]
          if (any(cond_vars %in% latent)) return(FALSE)
        }
        # Also skip if a latent is the focal pair variable (test[1] or test[2])
        # — these are internal MAG nodes and should not appear as observed test pairs
        if (test[1] %in% latent || test[2] %in% latent) return(FALSE)
        return(TRUE)
      },
      basis
    )
    n_removed <- n_before - length(basis)
    if (!quiet && n_removed > 0) {
      message(sprintf(
        "Removed %d untestable m-separation claim(s) where latent variable(s) appear in the conditioning set (Shipley & Douma 2021).",
        n_removed
      ))
    }
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
    hierarchical_info = hierarchical_info
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
        cat(format_dsep_test(test), "\n")
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
format_dsep_test <- function(test) {
  # Extract variables from formula
  # Use regex to strip out random terms like (1 | Group) before extracting vars
  test_str <- deparse(test)
  fixed_str <- gsub(
    "\\s*\\+\\s*\\(.*?\\|.*?\\)",
    "",
    paste(test_str, collapse = " ")
  )

  test_fixed <- stats::as.formula(fixed_str)
  vars <- all.vars(test_fixed)

  response <- as.character(test_fixed)[2]
  predictors <- setdiff(vars, response)

  if (length(predictors) == 0) {
    return(paste0("I( ", response, " , INVALID |  )"))
  }

  test_var <- predictors[1]

  if (length(predictors) == 1) {
    return(paste0("I( ", response, " , ", test_var, " |  )"))
  } else {
    # Sort for canonical representation
    cond_set <- paste(sort(predictors[-1]), collapse = ", ")
    return(paste0("I( ", response, " , ", test_var, " | ", cond_set, " )"))
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
plot_dsep.because <- function(object, ...) {
  if (is.null(object$dsep) || !object$dsep) {
    stop("plot_dsep requires a 'because' object fitted with dsep = TRUE.")
  }

  # summary.because handles the complex renaming and dummy variable matching
  s <- summary(object)

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
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "d-separation Independence Tests",
      subtitle = "Caterpillar plot of path coefficients with 95% Bayesian Credibility Intervals",
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
