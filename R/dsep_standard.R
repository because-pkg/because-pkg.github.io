# Standard DAG D-Separation
#
# Implements the Shipley (2000) basis-set approach using dagitty for models
# without latent variables.  Deterministic nodes (interactions, I() transforms)
# are retained as explicit graph nodes following Geiger, Verma & Pearl (1990).
#
# Primary function: dsep_standard()
#
# @keywords internal

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
        # Only inject if explicitly requested in equations OR formally declared as a structure
        is_requested <- any(sapply(random_terms, function(rt) rt$group == ir$group))
        is_structure <- !is.null(hierarchical_info$structure_levels) && 
                        (ir$group %in% names(hierarchical_info$structure_levels))
        
        if (is_requested || is_structure) {
          if (!any(sapply(combined_rand, function(x) x$group == ir$group))) {
            # Only add if the grouping variable is NOT the predictor itself
            test_var <- attr(t_eq, "test_var")
            if (is.null(test_var) || ir$group != test_var) {
              combined_rand[[length(combined_rand) + 1]] <- ir
            }
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

