#' Convert becauseR equations to ggm DAG adjacency matrix
#'
#' When \code{deterministic_terms} is supplied (a list returned by
#' \code{extract_deterministic_terms}), interaction and \code{I()} terms are
#' kept as **explicit intermediate nodes** in the DAG rather than collapsed to
#' their component variables.  This is required to produce the correct
#' conditional independence basis set following Geiger, Verma & Pearl (1990),
#' which extends d-separation to handle deterministic nodes.
#'
#' @param equations List of formulas
#' @param exclude_vars Character vector of variable names to exclude (e.g.,
#'   grouping variables)
#' @param deterministic_terms Optional named list returned by
#'   \code{extract_deterministic_terms}.  Each element must have
#'   \code{$original} (the R term string, e.g. \code{"BM:M"}) and
#'   \code{$internal_name} (the JAGS-safe node name, e.g. \code{"BM_x_M"}).
#' @return Named adjacency matrix in ggm format
#' @references
#'   Geiger, D., Verma, T., & Pearl, J. (1990). Identifying independence in
#'   Bayesian Networks. \emph{Networks}, 20(5), 507–534.
#' @keywords internal
equations_to_dag <- function(
    equations,
    exclude_vars = NULL,
    deterministic_terms = NULL
) {
    # Helper function to check if a term is a random effect term
    is_random_term <- function(term) {
        # Matches patterns like "(1|var)", "(1 | var)", etc.
        grepl("\\(.*\\|.*\\)", term) || grepl("^\\d+\\s*\\|\\s*\\w+$", term)
    }

    # Build a lookup: original R term string -> internal node name
    # e.g. "BM:M" -> "BM_x_M"
    det_lookup <- NULL
    if (!is.null(deterministic_terms) && length(deterministic_terms) > 0) {
        det_orig <- sapply(deterministic_terms, function(x) x$original)
        det_iname <- sapply(deterministic_terms, function(x) x$internal_name)
        det_lookup <- setNames(det_iname, det_orig)
    }

    # -----------------------------------------------------------------------
    # Collect all variable names, including deterministic node internal names
    # -----------------------------------------------------------------------
    all_vars <- unique(c(
        # Response-side variables
        sapply(equations, function(eq) as.character(eq[[2]])),
        # Predictor-side variables
        unlist(lapply(equations, function(eq) {
            term_labels <- attr(stats::terms(eq), "term.labels")
            vars <- character(0)
            for (term in term_labels) {
                if (is_random_term(term)) {
                    next
                }
                if (!is.null(det_lookup) && term %in% names(det_lookup)) {
                    # Deterministic term: add internal node name + its components
                    vars <- c(vars, det_lookup[[term]])
                    vars <- c(
                        vars,
                        all.vars(stats::as.formula(paste("~", term)))
                    )
                } else {
                    vars <- c(
                        vars,
                        all.vars(stats::as.formula(paste("~", term)))
                    )
                }
            }
            return(vars)
        }))
    ))

    # Exclude specified variables (e.g., grouping variables)
    if (!is.null(exclude_vars)) {
        all_vars <- setdiff(all_vars, exclude_vars)
    }

    n <- length(all_vars)
    dag <- matrix(0, n, n, dimnames = list(all_vars, all_vars))

    # -----------------------------------------------------------------------
    # Fill adjacency matrix
    # -----------------------------------------------------------------------
    for (eq in equations) {
        child <- as.character(eq[[2]])
        if (!child %in% all_vars) {
            next
        }

        term_labels <- attr(stats::terms(eq), "term.labels")

        # Identify random-effect term positions
        if (length(term_labels) > 0) {
            random_indices <- which(vapply(
                term_labels,
                is_random_term,
                FUN.VALUE = logical(1)
            ))
        } else {
            random_indices <- integer(0)
        }
        fixed_terms <- if (length(random_indices) > 0) {
            term_labels[-random_indices]
        } else {
            term_labels
        }

        for (term in fixed_terms) {
            if (!is.null(det_lookup) && term %in% names(det_lookup)) {
                # ---- Deterministic/interaction term ----
                # Route through the explicit deterministic node:
                #   component_var -> det_node -> child
                det_name <- det_lookup[[term]]

                if (det_name %in% all_vars) {
                    # Edge: det_node -> child
                    dag[det_name, child] <- 1

                    # Edges: each component variable -> det_node
                    comp_vars <- setdiff(
                        all.vars(stats::as.formula(paste("~", term))),
                        exclude_vars
                    )
                    for (cv in comp_vars) {
                        if (cv %in% all_vars) {
                            dag[cv, det_name] <- 1
                        }
                    }
                } else {
                    # det_name was excluded: fall back to direct edges
                    comp_vars <- setdiff(
                        all.vars(stats::as.formula(paste("~", term))),
                        exclude_vars
                    )
                    for (cv in comp_vars) {
                        if (cv %in% all_vars) dag[cv, child] <- 1
                    }
                }
            } else {
                # ---- Regular term ----
                term_vars <- tryCatch(
                    setdiff(
                        all.vars(stats::as.formula(paste("~", term))),
                        exclude_vars
                    ),
                    error = function(e) character(0)
                )
                for (v in term_vars) {
                    if (v %in% all_vars) dag[v, child] <- 1
                }
            }
        }
    }

    return(dag)
}

#' Extract bidirected edges (induced correlations) from MAG
#'
#' @param mag MAG adjacency matrix from DAG.to.MAG
#' @return List of variable pairs with induced correlations
#' @keywords internal
extract_bidirected_edges <- function(mag) {
    # Bidirected edges in MAG are coded as 100
    n <- nrow(mag)
    var_names <- rownames(mag)

    correlations <- list()

    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            if (mag[i, j] == 100 && mag[j, i] == 100) {
                # Bidirected edge found
                correlations[[length(correlations) + 1]] <- c(
                    var_names[i],
                    var_names[j]
                )
            }
        }
    }

    return(correlations)
}

#' Convert MAG basis set to becauseR formula format
#'
#' @param basis_set Basis set from basiSet.mag() or ggm::basiSet()
#' @param latent_children Optional character vector of variables that are direct children of latents
#' @param categorical_vars Optional named list of categorical variable info
#' @param family Optional named character vector of family distributions
#' @param deterministic_terms Optional named list from \code{extract_deterministic_terms}.
#'   Internal node names (e.g. \code{"BM_x_M"}) are rendered back to their
#'   original R syntax (e.g. \code{"BM:M"}) in the output formulas, following
#'   Geiger, Verma & Pearl (1990).
#' @param root_vars Optional character vector of root (exogenous) node names
#'   — nodes with no parents in the DAG.  If var1 is a root node and var2 is
#'   not, they are swapped so the non-root (downstream) variable becomes the
#'   response.  This enforces the convention that parent-only variables should
#'   always be predictors rather than responses in independence tests.
#' @return List of formulas with test_var attribute
#' @keywords internal
mag_basis_to_formulas <- function(
    basis_set,
    latent_children = NULL,
    categorical_vars = NULL,
    family = NULL,
    deterministic_terms = NULL,
    root_vars = NULL
) {
    # Build reverse lookup: internal_name -> original R syntax
    # e.g. "BM_x_M" -> "BM:M"
    det_rev_lookup <- NULL
    if (!is.null(deterministic_terms) && length(deterministic_terms) > 0) {
        det_iname <- sapply(deterministic_terms, function(x) x$internal_name)
        det_orig <- sapply(deterministic_terms, function(x) x$original)
        det_rev_lookup <- setNames(det_orig, det_iname)
    }

    # Helper: translate an internal node name back to displayable R syntax
    resolve_var <- function(v) {
        if (!is.null(det_rev_lookup) && v %in% names(det_rev_lookup)) {
            return(det_rev_lookup[[v]])
        }
        return(v)
    }

    if (is.null(basis_set) || length(basis_set) == 0) {
        return(list())
    }

    formulas <- list()

    for (test in basis_set) {
        # test is a vector: c(var1, var2, conditioning_vars...)
        var1 <- test[1]
        var2 <- test[2]
        cond_vars <- if (length(test) > 2) test[3:length(test)] else NULL

        # Apply ordering rule: if var1 is a latent child and var2 is not,
        # swap them so the latent child becomes the predictor (test variable).
        if (!is.null(latent_children)) {
            var1_is_latent_child <- var1 %in% latent_children
            var2_is_latent_child <- var2 %in% latent_children

            if (var1_is_latent_child && !var2_is_latent_child) {
                temp <- var1
                var1 <- var2
                var2 <- temp
            }
        }

        # Root-node rule: exogenous variables (no parents in the DAG) should
        # always be predictors, never responses.  If var1 is a root node and
        # var2 is not, swap so the non-root (downstream) variable is the response.
        if (!is.null(root_vars)) {
            var1_is_root <- var1 %in% root_vars
            var2_is_root <- var2 %in% root_vars

            if (var1_is_root && !var2_is_root) {
                temp <- var1
                var1 <- var2
                var2 <- temp
            }
        }

        # New Rule: Favor categorical variables as PREDICTORS (RHS/test_var)
        # If var1 (potential Response) is categorical, and var2 (Predictor) is NOT,
        # SWAP them so var1 becomes the predictor.
        # This allows us to use proper dummy variables in the test instead of
        # modeling the categorical variable as a Gaussian response.
        if (!is.null(categorical_vars)) {
            var1_is_cat <- var1 %in% names(categorical_vars)
            var2_is_cat <- var2 %in% names(categorical_vars)

            if (var1_is_cat && !var2_is_cat) {
                temp <- var1
                var1 <- var2
                var2 <- temp
            }
        }

        # FINAL RULE (User Override): Enforce directionality for species parameters.
        # Species parameters (p_, psi_, z_) should always be the RESPONSE in d-sep tests
        # with covariates, as covariates cannot cause parameters (biologically).
        # Regex to identify parameters: starts with p_, psi_, or z_
        is_param <- function(v) grepl("^(p_|psi_|z_)", v)

        if (is_param(var2) && !is_param(var1)) {
            # Predictor is parameter, Response is not. SWAP to make Parameter the Response.
            temp <- var1
            var1 <- var2
            var2 <- temp
        } else if (is_param(var2) && is_param(var1)) {
            # Both are parameters. Tie-breaker rule (User Request):
            # p_Species should be RESPONSE if paired with psi_Species/z_Species
            # Hierarchy: p_ > psi_/z_ > covariate

            is_p1 <- grepl("^p_", var1)
            is_p2 <- grepl("^p_", var2)

            # If Predictor is p_ and Response is NOT p_ (i.e. psi/z), SWAP.
            if (is_p2 && !is_p1) {
                temp <- var1
                var1 <- var2
                var2 <- temp
            }
        }

        # Build formula string

        # PREDICTOR RENAMING (User Request): Use psi_Species instead of z_Species/Species as predictor
        # This improves convergence for d-separation tests.
        if (!is.null(family)) {
            # Helper: rename if occupancy
            rename_if_occ <- function(v) {
                # Check if this variable is an occupancy variable
                # (either the base name or maybe already prefixed?)
                # Usually basis set uses base variable names (e.g. "Dingo")

                # Check if 'v' itself is in family as occupancy
                if (!is.na(family[v]) && family[v] == "occupancy") {
                    return(paste0("psi_", v))
                }

                # If it's already p_ or psi_ or z_, keep as is
                if (grepl("^(p_|psi_|z_)", v)) {
                    return(v)
                }

                return(v)
            }

            var2 <- rename_if_occ(var2)
            if (!is.null(cond_vars) && length(cond_vars) > 0) {
                cond_vars <- sapply(cond_vars, rename_if_occ)
            }
        }

        if (is.null(cond_vars) || length(cond_vars) == 0) {
            # Resolve internal names to R syntax for both var2 and cond_vars
            var2_r <- resolve_var(var2)
            formula_str <- paste(var1, "~", var2_r)
        } else {
            # Sort conditioning variables for a canonical representation and stable deduplication
            sorted_cond <- sort(cond_vars)
            # Resolve internal node names to original R syntax (e.g. BM_x_M -> BM:M)
            sorted_cond_r <- sapply(sorted_cond, resolve_var)
            var2_r <- resolve_var(var2)
            formula_str <- paste(
                var1,
                "~",
                var2_r,
                "+",
                paste(sorted_cond_r, collapse = " + ")
            )
        }
        f <- stats::as.formula(formula_str)
        attr(f, "test_var") <- var2_r # The variable being tested for independence

        formulas[[length(formulas) + 1]] <- f
    }

    # EXCLUSION RULE (User Request): Remove tests between p_Species and psi_Species (or Species itself)
    # for the same species. These are structurally coupled and testing independence is confusing/invalid.

    if (length(formulas) > 0) {
        keep_indices <- rep(TRUE, length(formulas))

        for (i in seq_along(formulas)) {
            f <- formulas[[i]]
            test_var_r <- attr(f, "test_var") # Original R syntax or internal name
            resp_var <- as.character(f)[2]

            resp_var <- as.character(f)[2]

            # -------------------------------------------------------------
            # 1. CYCLE PREVENTION: Detect Deterministic Descendants
            # -------------------------------------------------------------
            # If a conditioning variable is a deterministic term (e.g. BM:M)
            # and it contains the response or test_var as a component,
            # conditioning on it in a regression creates a JAGS cycle.
            if (!is.null(deterministic_terms)) {
                # [FIX] Use labels(terms(f)) instead of all.vars(f) to get the ACTUAL terms
                # on the RHS (including interactions and I(...) calls), otherwise they
                # are decomposed into components and the match fails.
                rhs_terms <- labels(terms(f))

                for (cv in rhs_terms) {
                    # Check if cv matches an original R syntax OR an internal JAGS name
                    det_match <- Filter(
                        function(dt) {
                            dt$original == cv || dt$internal_name == cv
                        },
                        deterministic_terms
                    )

                    if (length(det_match) > 0) {
                        # Extract components of the deterministic node
                        orig_cv <- det_match[[1]]$original
                        components <- all.vars(stats::as.formula(paste(
                            "~",
                            orig_cv
                        )))

                        orig_cv <- det_match[[1]]$original
                        components <- all.vars(stats::as.formula(paste(
                            "~",
                            orig_cv
                        )))

                        # Robust check: identify if EITHER the response OR the test variable
                        # is a parent/component of this deterministic conditioning node.
                        bad_components <- unique(c(
                            resp_var,
                            test_var_r,
                            sub("^psi_", "", resp_var),
                            sub("^psi_", "", test_var_r),
                            sub("^p_", "", resp_var),
                            sub("^p_", "", test_var_r)
                        ))

                        if (any(bad_components %in% components)) {
                            keep_indices[i] <- FALSE
                            message(paste0(
                                "  (Safety) Dropping cyclic d-sep test: ",
                                resp_var,
                                " ~ ",
                                test_var_r,
                                " | ",
                                cv,
                                " (conditioned on own component)"
                            ))
                            break
                        }
                    }
                }
            }
            if (!keep_indices[i]) {
                next
            }

            # -------------------------------------------------------------
            # 2. OCCUPANCY PARAMETER COUPLING
            # -------------------------------------------------------------
            # Helper to extract species name from p_Species or psi_Species or Species
            get_species <- function(x) {
                x <- sub("^p_", "", x)
                x <- sub("^psi_", "", x)
                x <- sub("^z_", "", x)
                return(x)
            }

            s1 <- get_species(resp_var)
            s2 <- get_species(test_var_r)

            # Check if one is p_ and the other is a state variable (psi, z, or raw species name)
            is_p1 <- grepl("^p_", resp_var)
            is_p2 <- grepl("^p_", test_var_r)

            is_state1 <- !is_p1
            is_state2 <- !is_p2

            # Condition: Same species AND one is p, one is state
            if (s1 == s2 && ((is_p1 && is_state2) || (is_state1 && is_p2))) {
                keep_indices[i] <- FALSE
            }
        }

        formulas <- formulas[keep_indices]
    }

    return(formulas)
}
