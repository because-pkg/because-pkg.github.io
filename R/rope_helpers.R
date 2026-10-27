#' Default and Resolved Region of Practical Equivalence (ROPE)
#'
#' Internal helpers for resolving family- and variable-specific ROPE intervals,
#' formatting human-readable explanations, and providing principled defaults.
#'
#' @keywords internal

#' Get default ROPE for a distribution family
#' @param fam Character; family name (e.g. "gaussian", "poisson", "binomial")
#' @return Numeric vector of length 2: c(min, max)
#' @keywords internal
get_family_default_rope <- function(fam) {
    if (is.null(fam) || is.na(fam)) fam <- "gaussian"
    if (inherits(fam, "because_family")) fam <- fam$name
    fam <- tolower(as.character(fam))

    if (fam %in% c("gaussian", "normal")) {
        # Kruschke standard: ±0.1 SD for standardized continuous traits
        c(-0.10, 0.10)
    } else if (fam %in% c("poisson", "negbinomial", "zip", "zinb")) {
        # ±5% relative rate ratio [0.95, 1.05]
        c(-0.05, 0.05)
    } else if (fam %in% c("binomial", "bernoulli", "ordinal", "multinomial")) {
        # ±0.1 latent SD on logit scale (0.1 * pi / sqrt(3) ≈ 0.181, OR [0.83, 1.20])
        c(-0.18, 0.18)
    } else {
        c(-0.10, 0.10)
    }
}

#' Resolve ROPE interval for a specific parameter
#'
#' Resolves in priority:
#' 1. Exact parameter match (e.g. rope[["beta_abundance_prey"]])
#' 2. Response variable match (e.g. rope[["abundance"]])
#' 3. Distribution family match (e.g. rope[["poisson"]])
#' 4. List default (rope[["default"]])
#' 5. Principled family default (get_family_default_rope)
#'
#' @param pname Parameter name
#' @param resp Response variable name
#' @param fam Family name
#' @param rope_spec User-supplied ROPE specification (NULL, TRUE, vector, or list)
#' @return Numeric vector of length 2
#' @keywords internal
resolve_rope_for_param <- function(pname, resp = NULL, fam = NULL, rope_spec = NULL) {
    # If rope_spec is a single numeric vector of length 2, use globally
    if (is.numeric(rope_spec) && length(rope_spec) == 2) {
        return(as.numeric(rope_spec))
    }

    # If rope_spec is a list or named vector, check hierarchy
    if (!is.null(rope_spec) && (is.list(rope_spec) || (!is.null(names(rope_spec)) && length(rope_spec) > 2))) {
        # 1. Exact parameter match
        if (!is.null(pname) && pname %in% names(rope_spec) && length(rope_spec[[pname]]) == 2) {
            return(as.numeric(rope_spec[[pname]]))
        }
        # 2. Response variable match
        if (!is.null(resp) && resp %in% names(rope_spec) && length(rope_spec[[resp]]) == 2) {
            return(as.numeric(rope_spec[[resp]]))
        }
        # 3. Family match
        fam_clean <- if (!is.null(fam)) tolower(as.character(fam)) else NULL
        if (!is.null(fam_clean) && fam_clean %in% names(rope_spec) && length(rope_spec[[fam_clean]]) == 2) {
            return(as.numeric(rope_spec[[fam_clean]]))
        }
        # 4. Default key in list
        if ("default" %in% names(rope_spec) && length(rope_spec[["default"]]) == 2) {
            return(as.numeric(rope_spec[["default"]]))
        }
    }

    # Fall back to family-specific default
    get_family_default_rope(fam)
}

#' Format human-readable explanations of ROPE intervals
#'
#' @param resolved_df Data frame containing Response, Family, r_min, r_max
#' @param data Optional model data frame to check if Gaussian traits are standardized
#' @return List with 'lines' (character vector) and 'note' (optional string)
#' @keywords internal
explain_rope_spec <- function(resolved_df, data = NULL) {
    if (is.null(resolved_df) || nrow(resolved_df) == 0) return(list(lines = character(0), note = NULL))

    # Keep all rows with non-NA values
    resolved_df <- unique(resolved_df[, c("Response", "Family", "r_min", "r_max")])

    # Find unique specs by Family, r_min, r_max
    u_specs <- unique(resolved_df[, c("Family", "r_min", "r_max")])
    total_responses <- unique(resolved_df$Response)
    n_total_resp <- length(total_responses)

    lines <- character(0)
    unstandardized_vars <- character(0)

    for (i in seq_len(nrow(u_specs))) {
        fam  <- tolower(as.character(u_specs$Family[i]))
        rmin <- u_specs$r_min[i]
        rmax <- u_specs$r_max[i]

        matching_rows <- resolved_df[tolower(as.character(resolved_df$Family)) == fam &
                                     resolved_df$r_min == rmin &
                                     resolved_df$r_max == rmax, ]
        resps <- unique(matching_rows$Response)
        resps <- resps[!is.na(resps) & resps != ""]

        fam_title <- if (fam %in% c("gaussian", "normal")) "Gaussian"
                     else if (fam %in% c("poisson", "negbinomial", "zip", "zinb")) {
                         if (fam == "poisson") "Poisson" else toupper(fam)
                     } else {
                         paste0(toupper(substr(fam, 1, 1)), substr(fam, 2, nchar(fam)))
                     }

        # Build response qualifier
        # If all responses in the plot share this exact same spec and there are multiple responses,
        # or if resp == "All", omit per-response listing
        resp_qualifier <- ""
        if ((nrow(u_specs) == 1 && n_total_resp > 1) || (length(resps) == 1 && resps[1] == "All")) {
            resp_qualifier <- ""
        } else if (length(resps) == 1) {
            resp_qualifier <- sprintf(" (%s)", resps[1])
        } else if (length(resps) <= 3) {
            resp_qualifier <- sprintf(" (%s)", paste(resps, collapse = ", "))
        } else {
            resp_qualifier <- sprintf(" (%s, +%d more)", paste(resps[1:2], collapse = ", "), length(resps) - 2)
        }

        if (fam %in% c("gaussian", "normal")) {
            lines <- c(lines, sprintf("%s%s: [%.2f, %.2f] (±%.2f SD if standardized)",
                                      fam_title, resp_qualifier, rmin, rmax, max(abs(c(rmin, rmax)))))
            # Check if trait is standardized in data
            if (!is.null(data)) {
                for (r in resps) {
                    if (r %in% names(data) && is.numeric(data[[r]])) {
                        s_val <- stats::sd(data[[r]], na.rm = TRUE)
                        if (!is.na(s_val) && (s_val > 1.5 || s_val < 0.6)) {
                            unstandardized_vars <- c(unstandardized_vars, sprintf("'%s' (SD = %.2f)", r, s_val))
                        }
                    }
                }
            }
        } else if (fam %in% c("poisson", "negbinomial", "zip", "zinb")) {
            irr_low <- exp(rmin)
            irr_upp <- exp(rmax)
            pct_low <- (irr_low - 1) * 100
            pct_upp <- (irr_upp - 1) * 100
            lines <- c(lines, sprintf("%s%s: [%.2f, %.2f] (Rate Ratio [%.2f, %.2f], %.1f%% to +%.1f%% relative count)",
                                      fam_title, resp_qualifier, rmin, rmax, irr_low, irr_upp, pct_low, pct_upp))
        } else if (fam %in% c("binomial", "bernoulli", "ordinal", "multinomial")) {
            or_low <- exp(rmin)
            or_upp <- exp(rmax)
            latent_sd <- max(abs(c(rmin, rmax))) / (pi / sqrt(3))
            lines <- c(lines, sprintf("%s%s: [%.2f, %.2f] (Odds Ratio [%.2f, %.2f], ±%.2f latent SD)",
                                      fam_title, resp_qualifier, rmin, rmax, or_low, or_upp, latent_sd))
        } else {
            lines <- c(lines, sprintf("%s%s: [%.2f, %.2f]", fam_title, resp_qualifier, rmin, rmax))
        }
    }

    note <- NULL
    if (length(unstandardized_vars) > 0) {
        note <- sprintf("Note: Variable(s) %s are not standardized. ROPE is evaluated in raw units; consider standardizing continuous variables or supplying a custom ROPE interval.",
                        paste(unique(unstandardized_vars), collapse = ", "))
    }

    list(lines = lines, note = note)
}
