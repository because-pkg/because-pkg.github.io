# DAG String Builder
#
# Converts a list of R model formulas to a dagitty-compatible DAG string,
# following the IDAG convention (Attia, Holliday & Oldmeadow 2022) where
# interaction terms (X:Y) and I() polynomial terms are explicit diamond nodes.
#
# Primary function: equations_to_dag_string()
#
# @keywords internal

#' Convert Equations List to DAGitty String
#'
#' Implements the Attia, Holliday & Oldmeadow (2022) IDAG convention:
#' interaction terms (e.g. \code{BM:M}) and \code{I()} terms are represented as
#' explicit intermediate nodes rather than collapsed to direct component->response
#' edges.
#'
#' @param equations List of formulas
#' @param induced_cors List of character vectors (pairs) for bidirected edges
#' @param family Optional named character vector of family distributions
#' @return A named list:
#'   \itemize{
#'     \item \code{dag_string} — dagitty-compatible DAG string
#'     \item \code{interaction_nodes} — named list: internal_name -> display label
#'   }
#' @references
#'   Attia, J., Holliday, E., & Oldmeadow, C. (2022). A proposal for capturing
#'   interaction and effect modification using DAGs.
#'   \emph{International Journal of Epidemiology}, 51(4), 1047--1053.
#' @noRd
equations_to_dag_string <- function(
    equations,
    induced_cors = NULL,
    family = NULL,
    poly_terms = NULL, # list from get_all_polynomial_terms; used for fitted models
    # where equations are already expanded (I(age^2) -> age_pow2)
    collapse_expanded = FALSE,
    extra_edges = character()
) {
    edges <- c()
    interaction_nodes <- list() # internal_name -> display label (e.g. "BM\u00d7M")

    # Build lookup: internal_name -> poly term info, for reconstructing diamond
    # nodes when equations are already expanded (fitted model path).
    poly_lookup <- list()
    if (!is.null(poly_terms)) {
        for (pt in poly_terms) {
            poly_lookup[[pt$internal_name]] <- pt
        }
    }

    # Helper: make a dagitty-safe node name from a term string
    # make_internal_name must match sanitize_term_name() from deterministic_nodes.R
    # so that node names in the plot align with JAGS parameter names like beta_weight_g_age_pow2.
    make_internal_name <- function(term) sanitize_term_name(term)

    # Helper: human-readable display label
    # For interactions: BM:M   -> BM×M
    # For I() powers:   I(age^2) -> age²  I(x^3) -> x³
    # For other I():    I(x+y)  -> I(x+y)  (keep as-is)
    make_display_label <- function(term) {
        if (grepl("^I\\(", term)) {
            # Check for simple power pattern: I(var^N)
            m <- regmatches(
                term,
                regexpr("^I\\(([a-zA-Z_][a-zA-Z0-9_]*)\\^([0-9]+)\\)$", term)
            )
            if (length(m) > 0) {
                inner <- sub("^I\\((.*)\\)$", "\\1", term)
                base <- sub("\\^.*$", "", inner)
                exp_n <- sub("^.*\\^", "", inner)
                superscripts <- c(
                    "\u2070",
                    "\u00b9",
                    "\u00b2",
                    "\u00b3",
                    "\u2074",
                    "\u2075",
                    "\u2076",
                    "\u2077",
                    "\u2078",
                    "\u2079"
                )
                n <- as.integer(exp_n)
                if (!is.na(n) && n >= 0 && n <= 9) {
                    return(paste0(base, superscripts[n + 1]))
                }
            }
            return(term) # fallback: keep as-is for complex I() expressions
        }
        gsub(":", "\u00d7", term) # "BM:M" -> "BM\u00d7M"
    }

    for (eq in equations) {
        resp <- all.vars(eq)[1]
        trm_lbls <- attr(terms(eq), "term.labels")

        # Optional: Collapse expanded names back to base names
        if (collapse_expanded) {
           resp <- gsub("(_[A-Za-z0-9]+|\\[[0-9]+\\])$", "", resp)
        }

        actual_resp <- resp

        if (length(trm_lbls) == 0) {
            next
        } # intercept-only, nothing to draw

        # Detect pure deterministic declaration: entire RHS is a single I() call.
        # e.g.  AgeClass ~ I(0 * (age < 0.02) + 1 * (age >= 0.02))
        # Mark the LHS itself as the deterministic/interaction node and draw edges
        # from the component variables directly — avoids a giant intermediate node.
        is_pure_det <- length(trm_lbls) == 1 && grepl("^I\\(", trm_lbls[1])

        if (is_pure_det) {
            term <- trm_lbls[1]
            interaction_nodes[[actual_resp]] <- actual_resp # LHS is the det. node
            components <- all.vars(stats::as.formula(paste("~", term)))
            for (comp in components) {
                edges <- c(edges, paste(actual_resp, "<-", comp))
            }
            next
        }

        for (term in trm_lbls) {
            # Skip random effects terms (e.g. 1 | year) which dagitty can't parse
            if (grepl("|", term, fixed = TRUE)) next
            
            # Optional: Collapse expanded names back to base names
            clean_term <- term
            if (collapse_expanded) {
               clean_term <- gsub("(_[A-Za-z0-9]+|\\[[0-9]+\\])$", "", term)
            }
            
            # Avoid self-loops if resp and term matched the same base name
            if (clean_term == resp) next

            is_interaction <- grepl(":", clean_term, fixed = TRUE) &&
                !grepl("^I\\(", clean_term)
            is_I_call <- grepl("^I\\(", clean_term)

            if (is_interaction || is_I_call) {
                # --- Deterministic / interaction node (within a mixed equation) ---
                # X:Y interactions and I() polynomial terms get a diamond node.
                iname <- make_internal_name(clean_term)
                d_label <- make_display_label(clean_term)
                interaction_nodes[[iname]] <- d_label

                # Component variables: split on : for interactions, all.vars for I()
                if (is_interaction) {
                    components <- strsplit(clean_term, ":", fixed = TRUE)[[1]]
                } else {
                    components <- all.vars(stats::as.formula(paste("~", clean_term)))
                }

                for (comp in components) {
                    edges <- c(edges, paste(iname, "<-", comp))
                }
                # Route through the deterministic node: response <- iname
                # For X:Y interactions, components may not appear as separate
                # explicit predictors, so also add comp -> resp edges.
                # For I() polynomial terms (e.g. I(age^2)), the base variable
                # (age) is typically an explicit separate term in the same equation,
                # so comp -> resp is already added by the regular-predictor branch.
                # Adding it again here creates a redundant hidden edge that
                # collides with the direct age -> weight_g arrow in the layout.
                edges <- c(edges, paste(actual_resp, "<-", iname))
                if (is_interaction) {
                    for (comp in components) {
                        edges <- c(edges, paste(actual_resp, "<-", comp))
                    }
                }
            } else {
                # --- Regular predictor ---
                # If we are collapsing, we skip the expansion-specific diamond nodes
                # and just draw an edge from the base variable to the response.
                if (collapse_expanded && term %in% names(poly_lookup)) {
                   # Skip the diamond node, handled by 'clean_term' falling through to Standard path
                } else if (!collapse_expanded && term %in% names(poly_lookup)) {
                    pt <- poly_lookup[[term]]
                    iname <- pt$internal_name # e.g. age_pow2
                    # ... (rest of diamond logic) ...
                    superscripts <- c("\u2070","\u00b9","\u00b2","\u00b3","\u2074","\u2075","\u2076","\u2077","\u2078","\u2079")
                    n <- as.integer(pt$power)
                    d_label <- if (!is.na(n) && n >= 0 && n <= 9) paste0(pt$base_var, superscripts[n+1]) else iname
                    interaction_nodes[[iname]] <- d_label
                    edges <- c(edges, paste(iname, "<-", pt$base_var))
                    edges <- c(edges, paste(actual_resp, "<-", iname))
                } else {
                    edges <- c(edges, paste(actual_resp, "<-", clean_term))
                }
            }
        }
    }

    # Append any extra edges injected by extension packages via dag_expand_hook
    if (length(extra_edges) > 0) {
        edges <- c(edges, extra_edges)
    }

    if (!is.null(induced_cors)) {
        for (pair in induced_cors) {
            if (length(pair) == 2) {
                edges <- c(edges, paste(pair[1], "<->", pair[2]))
            }
        }
    }

    list(
        dag_string = paste("dag {", paste(unique(edges), collapse = "; "), "}"),
        interaction_nodes = interaction_nodes
    )
}
