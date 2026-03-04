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
#' @param quiet Logical; if FALSE (default), print the basis set and MAG structure.
#'   If TRUE, suppress informational output.
#' @param random_terms Optional list of random effects (group, type) parsed from equations.
#' @param hierarchical_info Internal argument used to pass data hierarchy information
#'   (levels, grouping variables) for future implementation of multilevel d-separation
#'   tests (following Shipley 2009). Currently unused by the d-separation logic.
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

# Standard d-separation for DAGs (using ggm::basiSet)
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

  # Use ggm::basiSet to get the correct d-separation basis set
  if (!requireNamespace("ggm", quietly = TRUE)) {
    stop("Package 'ggm' is required for d-separation tests.")
  }

  basis <- ggm::basiSet(dag)

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

      # Find random terms for response (handle both Species and psi_Species)
      base_resp <- sub("^psi_", "", resp)
      vocab_rand <- Filter(
        function(x) x$response == resp || x$response == base_resp,
        random_terms
      )

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

  # Always suppress DAG.to.MAG output
  invisible(capture.output(
    {
      mag <- suppressMessages(DAG.to.MAG(dag, latents = latent))
    },
    type = "output"
  ))

  # Extract basis set from MAG
  invisible(capture.output(
    {
      basis <- suppressMessages(basiSet.mag(mag))
    },
    type = "output"
  ))

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
        function(x) x$response == resp || x$response == base_resp,
        random_terms
      )

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
  correlations <- extract_bidirected_edges(mag)

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
    mag = mag
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
