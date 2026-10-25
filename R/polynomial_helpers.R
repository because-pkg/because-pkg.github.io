#' Extract Polynomial Terms from Formula
#'
#' Detects I(variable^power) patterns in formula and extracts transformation info
#'
#' @param formula Formula object
#' @return List with polynomial term information, or NULL if none found
#' @keywords internal
extract_polynomial_terms <- function(formula) {
    # Convert formula to character
    formula_str <- paste(deparse(formula), collapse = " ")

    # Pattern to match I(var^power)
    # Matches: I(age^2), I(weight^3), etc.
    poly_pattern <- "I\\(\\s*([a-zA-Z_][a-zA-Z0-9_]*)\\s*\\^\\s*([0-9]+)\\s*\\)"

    # Find all matches
    matches <- gregexpr(poly_pattern, formula_str, perl = TRUE)

    if (matches[[1]][1] == -1) {
        return(NULL) # No polynomial terms found
    }

    # Extract matched strings
    match_strings <- regmatches(formula_str, matches)[[1]]

    if (length(match_strings) == 0) {
        return(NULL)
    }

    # Parse each polynomial term
    poly_terms <- list()

    for (match_str in match_strings) {
        # Extract variable and power
        parts <- regmatches(
            match_str,
            regexec(poly_pattern, match_str, perl = TRUE)
        )[[1]]

        if (length(parts) >= 3) {
            base_var <- parts[2]
            power <- as.integer(parts[3])

            # Generate internal name
            internal_name <- paste0(base_var, "_pow", power)

            poly_terms[[length(poly_terms) + 1]] <- list(
                original = match_str,
                base_var = base_var,
                power = power,
                internal_name = internal_name
            )
        }
    }

    if (length(poly_terms) == 0) {
        return(NULL)
    }

    return(poly_terms)
}


#' Expand Polynomial Terms in Formula
#'
#' Replaces I(var^power) with internal variable names
#'
#' @param formula Formula object
#' @param poly_terms List of polynomial term info from extract_polynomial_terms
#' @return Modified formula with internal names
#' @keywords internal
expand_polynomial_formula <- function(formula, poly_terms) {
    if (is.null(poly_terms)) {
        return(formula)
    }

    formula_str <- paste(deparse(formula), collapse = " ")

    # Replace each I(var^power) with internal name
    for (term in poly_terms) {
        # Escape special regex characters in original
        original_escaped <- gsub("([\\^\\(\\)])", "\\\\\\1", term$original)
        formula_str <- gsub(original_escaped, term$internal_name, formula_str)
    }

    # Convert back to formula
    as.formula(formula_str)
}


#' Generate BUGS Code for Polynomial Transformations
#'
#' Creates deterministic nodes for polynomial terms in JAGS and NIMBLE models
#'
#' @param poly_terms List of polynomial term info
#' @return Character vector of BUGS code lines
#' @keywords internal
generate_polynomial_jags <- function(poly_terms) {
    if (is.null(poly_terms) || length(poly_terms) == 0) {
        return(character(0))
    }

    jags_lines <- character(length(poly_terms))

    for (i in seq_along(poly_terms)) {
        term <- poly_terms[[i]]

        # Generate: var_pow2[i] <- var[i]^2
        jags_lines[i] <- sprintf(
            "    %s[i] <- %s[i]^%d",
            term$internal_name,
            term$base_var,
            term$power
        )
    }

    jags_lines
}


#' Get All Polynomial Terms from Equations
#'
#' Extracts all unique polynomial terms across all equations
#'
#' @param equations List of formula objects
#' @return List of all polynomial terms
#' @keywords internal
get_all_polynomial_terms <- function(equations) {
    all_poly_terms <- list()

    for (eq in equations) {
        poly_terms <- extract_polynomial_terms(eq)

        if (!is.null(poly_terms)) {
            all_poly_terms <- c(all_poly_terms, poly_terms)
        }
    }

    # Remove duplicates (same base_var and power)
    if (length(all_poly_terms) == 0) {
        return(NULL)
    }

    # Create unique keys
    keys <- sapply(all_poly_terms, function(x) {
        paste(x$base_var, x$power, sep = "_")
    })

    unique_indices <- !duplicated(keys)
    all_poly_terms[unique_indices]
}
