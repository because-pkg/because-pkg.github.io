#' Extract Deterministic Terms from Formulas
#'
#' Scans a list of formulas for terms that require deterministic nodes in JAGS
#' (interactions, I() calls, etc.)
#'
#' @param equations List of formulas
#' @return List of deterministic term definitions
#' @keywords internal
extract_deterministic_terms <- function(equations) {
    # Pre-scan equations for deterministic assignments (e.g. X ~ I(A+B))
    # If a variable is defined purely by a deterministic term, we use the variable
    # name as the internal name for the node to avoid redundant intermediate nodes.
    assignments <- list()
    for (eq in equations) {
        resp <- as.character(formula(eq))[2]
        # Get RHS terms strictly (no random effects)
        rhs_terms <- attr(terms(formula(eq)), "term.labels")
        # Check for random effects in the formula itself to be safe
        has_random <- grepl("\\|", as.character(formula(eq))[3])
        
        if (!has_random && length(rhs_terms) == 1 && (grepl(":", rhs_terms) || grepl("\\(", rhs_terms))) {
            assignments[[rhs_terms]] <- resp
        }
    }

    terms_list <- list()

    for (eq in equations) {
        # Get all terms including interactions
        eq_terms <- attr(terms(formula(eq)), "term.labels")

        for (term in eq_terms) {
            # Skip random effects
            if (grepl("\\|", term)) next

            # Check if term needs deterministic handling
            if (
                grepl(":", term) || (grepl("\\(", term) && grepl("\\)", term))
            ) {
                # Use assigned variable name if this term is a pure assignment
                if (term %in% names(assignments)) {
                    internal_name <- assignments[[term]]
                } else {
                    internal_name <- sanitize_term_name(term)
                }

                if (!internal_name %in% names(terms_list)) {
                    terms_list[[internal_name]] <- list(
                        original = term,
                        internal_name = internal_name,
                        expression = term_to_jags_expression(term)
                    )
                }
            }
        }
    }

    return(terms_list)
}

#' Sanitize Term Name for JAGS
#'
#' Converts complex R terms into valid JAGS variable names
#'
#' @param term Character string (e.g., "A:B", "I(A^2)")
#' @return Sanitized string (e.g., "A_x_B", "A_pow2")
#' @keywords internal
sanitize_term_name <- function(term) {
    # 1. Interactions A:B -> A_x_B
    out <- gsub(":", "_x_", term)

    # 2. Remove I(...) wrapper
    if (grepl("^I\\(", out)) {
        out <- sub("^I\\((.*)\\)$", "\\1", out)
    }

    # 3. Handle powers ^ -> _pow
    out <- gsub("\\^", "_pow", out)

    # 4. Handle logical AND/OR
    # R uses & / | but JAGS internal naming needs to be clean
    out <- gsub("&", "_and_", out)
    out <- gsub("\\|", "_or_", out)

    # 5. Handle comparisons (simple replacement)
    out <- gsub(">", "_gt_", out)
    out <- gsub("<", "_lt_", out)
    out <- gsub("==", "_eq_", out)
    out <- gsub(">=", "_gte_", out)
    out <- gsub("<=", "_lte_", out)

    # 5. Handle arithmetic
    out <- gsub("\\+", "_plus_", out)
    # Be careful with minus vs hyphen in variable names
    out <- gsub("\\-", "_minus_", out)
    out <- gsub("\\*", "_times_", out)
    out <- gsub("/", "_div_", out)

    # 6. Clean up any remaining non-alphanumeric chars
    out <- gsub("[^a-zA-Z0-9_]", "_", out)

    # 7. Remove duplicate underscores
    out <- gsub("_+", "_", out)

    # 8. Remove leading/trailing underscores
    out <- gsub("^_", "", out)
    out <- gsub("_$", "", out)

    # Check if name is complex (arbitrary heuristic: > 20 chars or contains 'times'/'plus')
    if (nchar(out) > 32) {
        # Extract vars safely
        vars <- tryCatch(all.vars(parse(text = term)), error = function(e) {
            character(0)
        })

        if (length(vars) > 0) {
            # Create base name from vars (e.g. "Age" or "Age_Sex")
            base <- paste(head(vars, 3), collapse = "_") # Limit to first 3 vars

            # Create hash for uniqueness (sum of ascii codes)
            checksum <- sprintf("%x", sum(utf8ToInt(out)))

            # New readable name
            out <- paste0("det_", base, "_", checksum)
        } else {
            # Fallback if no vars found
            # Truncate existing if too long
            if (nchar(out) > 25) {
                checksum <- sprintf("%x", sum(utf8ToInt(out)))
                out <- paste0(substr(out, 1, 25), "_", checksum)
            }

            if (grepl("^[0-9]", out)) {
                out <- paste0("det_", out)
            }
        }
    } else {
        # Valid short name, ensure prefix
        if (grepl("^[0-9]", out)) {
            out <- paste0("det_", out)
        }
    }

    return(out)
}

#' Convert R Term to JAGS Expression
#'
#' Transforms R syntax into JAGS-compatible deterministic code
#'
#' @param term Original R term (e.g., "A:B")
#' @return JAGS expression (e.g., "A\\[i] * B\\[i]")
#' @keywords internal
term_to_jags_expression <- function(term) {
    # 1. Handle Interaction A:B -> A[i] * B[i]
    if (grepl(":", term) && !grepl("\\(", term)) {
        parts <- unlist(strsplit(term, ":"))
        # Add [i] to each part
        parts_indexed <- paste0(parts, "[i]")
        return(paste(parts_indexed, collapse = " * "))
    }

    # 2. Handle I(...) wrapper - strip it for parsing
    # Remove newlines for easier parsing
    term <- gsub("\\n", " ", term)
    term <- gsub("\\s+", " ", term) # Collapse multiple spaces

    # 2. Handle I(...) wrapper - strip it for parsing
    inner <- term
    if (grepl("^I\\(", term)) {
        # Use simple string extraction to be safer than regex
        inner <- substr(term, 3, nchar(term) - 1)
    }

    # Split by non-word chars
    tokens <- strsplit(inner, "([\\+\\-\\*\\/\\^\\<\\>\\=\\(\\)\\s])")[[1]]
    delimiters <- strsplit(inner, "([a-zA-Z0-9_\\.]+)")[[1]]

    # handle basic variables in formulas:
    expression <- inner

    # Get all variables in the expression using all.vars
    # Note: all.vars(parse(text="A*B")) -> c("A", "B")
    vars <- all.vars(parse(text = inner))

    # Sort by length (descending) to avoid partial replacement (e.g. replacing 'Age' inside 'Age2')
    vars <- vars[order(nchar(vars), decreasing = TRUE)]

    for (v in vars) {
        # Replace variable with variable[i]
        # Use simple word boundary
        pattern <- paste0("\\b", v, "\\b")
        replacement <- paste0(v, "[i]")
        expression <- gsub(pattern, replacement, expression)
    }

    if (grepl("==", expression)) {
        # Simple regex for A == B -> equals(A, B)
        expression <- gsub(
            "([^[:space:]()=]+)[[:space:]]*==[[:space:]]*([^[:space:]()=]+)",
            "equals(\\1, \\2)",
            expression
        )
    }

    # Handle logical operators
    # Replace runs of & with && and | with ||
    expression <- gsub("&+", "&&", expression)
    expression <- gsub("\\|+", "||", expression)

    return(expression)
}

#' Generate JAGS Code for Deterministic Nodes
#'
#' writers the "Internal_Var\\[i] <- ..." lines
#'
#' @param deterministic_terms List from extract_deterministic_terms
#' @return Vector of JAGS code strings
#' @keywords internal
generate_deterministic_jags <- function(deterministic_terms) {
    if (length(deterministic_terms) == 0) {
        return(character(0))
    }

    lines <- c()
    for (term in deterministic_terms) {
        line <- sprintf("    %s[i] <- %s", term$internal_name, term$expression)
        lines <- c(lines, line)
    }
    return(lines)
}
