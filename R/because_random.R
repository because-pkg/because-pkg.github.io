#' Extract Random Effects from Formulas
#'
#' Parses model formulas to separate fixed effects from random effects terms
#' like (1|Group).
#'
#' @param equations A list of formulas
#' @return A list containing:
#'   \item{fixed_equations}{List of formulas with random terms removed}
#'   \item{random_terms}{List of extracted random effect specs (group, type)}
#' @noRd
extract_random_effects <- function(equations) {
    fixed_equations <- list()
    random_terms <- list()

    for (i in seq_along(equations)) {
        eq <- equations[[i]]
        response <- All.vars(eq)[1]

        # Get terms from formula
        tf <- terms(eq)
        term_labels <- attr(tf, "term.labels")

        # Parse RHS to separate random effects from fixed effects using string manipulation
        # This avoids issues with terms() parsing non-standard (1|Group) syntax
        char_eq <- as.character(eq)
        lhs <- char_eq[2]
        rhs <- char_eq[3]

        rhs_parts <- trimws(strsplit(rhs, "\\+")[[1]])

        fixed_parts <- c()

        for (part in rhs_parts) {
            # Use helper to extract term
            term <- parse_random_part(part)

            if (!is.null(term)) {
                # Add to random terms
                random_terms[[length(random_terms) + 1]] <- list(
                    response = lhs,
                    group = term$group,
                    type = term$type
                )
            } else {
                fixed_parts <- c(fixed_parts, part)
            }
        }

        # Reconstruct fixed formula
        if (length(fixed_parts) == 0) {
            fixed_parts <- "1"
        }
        new_rhs <- paste(fixed_parts, collapse = " + ")
        fixed_equations[[i]] <- as.formula(
            paste(lhs, "~", new_rhs),
            env = environment(eq)
        )
    }

    return(list(
        fixed_equations = fixed_equations,
        random_terms = random_terms
    ))
}

#' Create Group Structures for JAGS
#'
#' Creates indices and precision matrices for random effect groups.
#'
#' @param data Data frame or list
#' @param random_terms List of random terms extracted from equations
#' @return List of structure definitions to be added to data
#' @noRd
create_group_structures <- function(data, random_terms) {
    # We need unique groups across all equations
    # If Y ~ (1|ID) and Z ~ (1|ID), it's the same ID structure

    all_groups <- unique(sapply(random_terms, function(x) x$group))
    structures <- list()
    data_updates <- list()

    for (grp in all_groups) {
        if (!grp %in% names(data)) {
            stop(
                "Grouping variable '",
                grp,
                "' used in random effects not found in data."
            )
        }

        vals <- data[[grp]]

        # Handle factor/character
        if (is.factor(vals)) {
            levels <- levels(vals)
            indices <- as.integer(vals)
        } else {
            # Treat as character/factor
            vals <- as.factor(vals)
            levels <- levels(vals)
            indices <- as.integer(vals)
        }

        n_groups <- length(levels)

        # Create simple identity precision
        # In JAGS: u[1:N] ~ dmnorm(zeros, Prec)
        # Prec = Identity

        prec_mat <- diag(n_groups)
        # Precision matrix for JAGS
        # dimnames(prec_mat) <- list(levels, levels)

        # Store definition
        # Naming convention: group_{grp} for indices
        # Prec_{grp} for precision

        index_name <- paste0("group_", grp)
        prec_name <- paste0("Prec_", grp)
        zeros_name <- paste0("zeros_", grp)

        data_updates[[index_name]] <- indices
        data_updates[[prec_name]] <- prec_mat
        data_updates[[zeros_name]] <- rep(0, n_groups)
        data_updates[[paste0("N_", grp)]] <- n_groups

        # We add this group name to the list of structures
        structures[[grp]] <- list(
            name = grp,
            n = n_groups
        )
    }

    return(list(
        structures = structures,
        data_updates = data_updates
    ))
}

# Helper to safely all vars from list of formulas (if needed locally)
All.vars <- function(x) all.vars(x)


#' Expand Nesting Syntax in Random Effects
#'
#' Converts lme4-style nesting (1|A/B) to (1|A) + (1|A:B)
#'
#' @param rhs Character string, RHS of random effects formula
#' @return Expanded character string
#' @keywords internal
expand_nesting_syntax <- function(rhs) {
    # Pattern to match (1|A/B) or (1|A/B/C) etc.
    nesting_pattern <- "\\(1\\s*\\|\\s*([^)]+)\\)"

    # Find all random effect terms
    matches <- gregexpr(nesting_pattern, rhs)
    match_strings <- regmatches(rhs, matches)[[1]]

    if (length(match_strings) == 0) {
        return(rhs) # No nesting found
    }

    # Process each matched term
    for (match_str in match_strings) {
        # Extract the grouping part (everything after |)
        group_part <- sub("\\(1\\s*\\|\\s*(.+)\\)", "\\1", match_str)

        # Check if it contains nesting (/)
        if (grepl("/", group_part)) {
            # Split by /
            levels <- trimws(strsplit(group_part, "/")[[1]])

            # Build expanded terms
            # (1|A/B) -> (1|A) + (1|A:B)
            # (1|A/B/C) -> (1|A) + (1|A:B) + (1|A:B:C)
            expanded_terms <- character(length(levels))

            for (i in seq_along(levels)) {
                if (i == 1) {
                    expanded_terms[i] <- paste0("(1|", levels[i], ")")
                } else {
                    # Create interaction: A:B or A:B:C
                    interaction <- paste(levels[1:i], collapse = ":")
                    expanded_terms[i] <- paste0("(1|", interaction, ")")
                }
            }

            # Join with +
            replacement <- paste(expanded_terms, collapse = " + ")

            # Replace in original string (use escaping for special regex chars)
            match_escaped <- gsub("([\\(\\)\\|\\:])", "\\\\\\1", match_str)
            rhs <- sub(match_escaped, replacement, rhs)
        }
    }

    return(rhs)
}

# Parse Random Effect Part
parse_random_part <- function(part) {
    # Check if part contains a pipe character (random effect indicator)
    if (grepl("\\|", part)) {
        # Remove parentheses and leading ~ if present
        inner <- part
        # Strip leading tilde if present
        if (grepl("^~", inner)) {
            inner <- sub("^~", "", inner)
        }
        # Strip outer parentheses
        if (grepl("^\\(.*\\)$", inner)) {
            inner <- sub("^\\((.*)\\)$", "\\1", inner)
        }

        split_term <- strsplit(inner, "\\|")[[1]]

        if (length(split_term) != 2) {
            warning("Invalid random effect syntax: ", part)
            return(NULL)
        }

        lhs_term <- trimws(split_term[1])
        group_var <- trimws(split_term[2])

        if (lhs_term != "1") {
            warning(
                "Random slopes (",
                lhs_term,
                "|",
                group_var,
                ") are not yet implemented. ",
                "Defaulting to random intercept (1|",
                group_var,
                ")."
            )
        }

        return(list(
            group = group_var,
            type = "intercept"
        ))
    }
    return(NULL)
}

#' Parse Global Random Argument
#'
#' @param random_arg Formula (e.g. ~ (1|Site)) or character string
#' @param equations List of equations (to identify LHS variables)
#' @param all_vars Optional character vector of all variables to map to (for d-sep)
#' @return List of random terms (response, group, type)
#' @noRd
parse_global_random <- function(random_arg, equations, all_vars = NULL) {
    if (is.null(random_arg)) {
        return(list())
    }

    if (inherits(random_arg, "formula")) {
        char_arg <- as.character(random_arg)
        # Use last part which contains the RHS usually
        rhs <- char_arg[length(char_arg)]
    } else {
        rhs <- as.character(random_arg)
    }

    # Convert nesting syntax (1|A/B) to (1|A) + (1|A:B) to match lme4 behavior
    rhs <- expand_nesting_syntax(rhs)

    parts <- trimws(strsplit(rhs, "\\+")[[1]])
    global_terms <- list()

    for (part in parts) {
        term <- parse_random_part(part)
        if (!is.null(term)) {
            global_terms[[length(global_terms) + 1]] <- term
        }
    }

    final_terms <- list()
    if (length(global_terms) > 0) {
        # Identify variables to map to
        if (!is.null(all_vars)) {
            responses <- unique(all_vars)
        } else {
            responses <- unique(sapply(equations, function(eq) {
                as.character(formula(eq))[2]
            }))
        }

        for (resp in responses) {
            for (term in global_terms) {
                final_terms[[length(final_terms) + 1]] <- list(
                    response = resp,
                    group = term$group,
                    type = term$type
                )
            }
        }
    }

    return(final_terms)
}
