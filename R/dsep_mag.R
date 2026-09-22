# MAG M-Separation (Models with Latent Variables)
#
# Implements DAG-to-MAG conversion and m-separation basis-set generation
# following Shipley & Douma (2021).  Handles induced correlations (bidirected
# edges) arising from shared unmeasured common causes.
#
# Primary functions:
#   - dsep_with_latents()   MAG m-separation entry point
#   - format_dsep_test()    format a test formula for printing
#
# @keywords internal

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

