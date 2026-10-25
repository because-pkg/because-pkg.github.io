#' Extract Deterministic Terms from Formulas
#'
#' Scans a list of formulas for terms that require deterministic node definitions
#' (interactions, \code{I()} expressions, etc.) across model building, d-separation,
#' and DAG construction.
#'
#' @param equations List of formulas
#' @return List of deterministic term definitions
#' @keywords internal
extract_deterministic_terms <- function(equations) {
    terms_list <- list()

    for (eq in equations) {
        eq_form <- if (is.list(eq) && "formula" %in% names(eq)) eq$formula else formula(eq)
        eq_terms <- attr(terms(eq_form), "term.labels")

        for (term in eq_terms) {
            if (grepl("\\|", term)) next

            if (
                grepl(":", term) || (grepl("\\(", term) && grepl("\\)", term))
            ) {
                internal_name <- sanitize_term_name(term)

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

#' Sanitize Formula Term Name
#'
#' Converts complex R formula terms into valid alphanumeric node/parameter names
#' for BUGS engines (JAGS/NIMBLE), parameter maps, d-separation tests, and DAG layouts.
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
    out <- gsub("&", "_and_", out)
    out <- gsub("\\|", "_or_", out)

    # 5. Handle comparisons
    out <- gsub(">", "_gt_", out)
    out <- gsub("<", "_lt_", out)
    out <- gsub("==", "_eq_", out)
    out <- gsub(">=", "_gte_", out)
    out <- gsub("<=", "_lte_", out)

    # 6. Handle arithmetic
    out <- gsub("\\+", "_plus_", out)
    out <- gsub("\\-", "_minus_", out)
    out <- gsub("\\*", "_times_", out)
    out <- gsub("/", "_div_", out)

    # 7. Clean up any remaining non-alphanumeric chars
    out <- gsub("[^a-zA-Z0-9_]", "_", out)

    # 8. Remove duplicate and leading/trailing underscores
    out <- gsub("_+", "_", out)
    out <- gsub("^_", "", out)
    out <- gsub("_$", "", out)

    # Shorten overly long names (> 32 chars) using a deterministic hash
    if (nchar(out) > 32) {
        vars <- tryCatch(all.vars(parse(text = term)), error = function(e) {
            character(0)
        })

        if (length(vars) > 0) {
            base <- paste(head(vars, 3), collapse = "_")
            checksum <- sprintf("%x", sum(utf8ToInt(out)))
            out <- paste0("det_", base, "_", checksum)
        } else {
            if (nchar(out) > 25) {
                checksum <- sprintf("%x", sum(utf8ToInt(out)))
                out <- paste0(substr(out, 1, 25), "_", checksum)
            }

            if (grepl("^[0-9]", out)) {
                out <- paste0("det_", out)
            }
        }
    } else if (grepl("^[0-9]", out)) {
        out <- paste0("det_", out)
    }

    return(out)
}

#' Convert R Formula Term to BUGS Expression
#'
#' Transforms R formula syntax into indexed BUGS deterministic expressions for JAGS and NIMBLE.
#'
#' @param term Original R term (e.g., "A:B")
#' @return BUGS expression (e.g., "A\\[i] * B\\[i]")
#' @keywords internal
term_to_jags_expression <- function(term) {
    # 1. Handle Interaction A:B -> A[i] * B[i]
    if (grepl(":", term) && !grepl("\\(", term)) {
        parts <- unlist(strsplit(term, ":"))
        parts_indexed <- paste0(parts, "[i]")
        return(paste(parts_indexed, collapse = " * "))
    }

    # 2. Strip I(...) wrapper and normalize whitespace
    term <- gsub("\\n", " ", term)
    term <- gsub("\\s+", " ", term)

    inner <- term
    if (grepl("^I\\(", term)) {
        inner <- substr(term, 3, nchar(term) - 1)
    }

    expression <- inner
    vars <- all.vars(parse(text = inner))

    # Sort by length (descending) to avoid partial replacement (e.g. replacing 'Age' inside 'Age2')
    vars <- vars[order(nchar(vars), decreasing = TRUE)]

    for (v in vars) {
        pattern <- paste0("\\b", v, "\\b")
        replacement <- paste0(v, "[i]")
        expression <- gsub(pattern, replacement, expression)
    }

    if (grepl("==", expression)) {
        expression <- gsub(
            "([^[:space:]()=]+)[[:space:]]*==[[:space:]]*([^[:space:]()=]+)",
            "equals(\\1, \\2)",
            expression
        )
    }

    # Replace runs of & with && and | with ||
    expression <- gsub("&+", "&&", expression)
    expression <- gsub("\\|+", "||", expression)

    return(expression)
}
