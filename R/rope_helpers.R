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
#' @param standardized Logical. If TRUE, forces Gaussian standard deviation benchmark c(-0.10, 0.10).
#' @return Numeric vector of length 2
#' @keywords internal
resolve_rope_for_param <- function(pname, resp = NULL, fam = NULL, rope_spec = NULL, standardized = FALSE) {
    if (isTRUE(standardized)) {
        fam <- "gaussian"
    }

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

#' Extract numeric posterior draws for a d-separation test parameter
#'
#' Checks `idx_hint` first (preserving 1-to-1 ordering with `dsep$results`),
#' then falls back to scanning `dsep_results` by parameter name.
#'
#' @param dsep_results List of per-test result entries containing `$samples`
#' @param param_name Parameter name to match
#' @param idx_hint Optional integer index hint
#' @return Numeric vector of posterior draws, or `NULL` if not found
#' @keywords internal
.find_dsep_param_draws <- function(dsep_results, param_name, idx_hint = NULL) {
    if (is.null(dsep_results) || is.null(param_name) || is.na(param_name)) {
        return(NULL)
    }

    extract_from_entry <- function(dr) {
        if (is.null(dr) || is.null(dr$samples)) return(NULL)
        samps <- as.matrix(dr$samples)
        if (param_name %in% colnames(samps)) {
            return(as.numeric(samps[, param_name]))
        }
        NULL
    }

    if (!is.null(idx_hint) && idx_hint >= 1L && idx_hint <= length(dsep_results)) {
        vec <- extract_from_entry(dsep_results[[idx_hint]])
        if (!is.null(vec)) return(vec)
    }

    for (dr in dsep_results) {
        vec <- extract_from_entry(dr)
        if (!is.null(vec)) return(vec)
    }

    NULL
}

#' Format concise response-name qualifier for ROPE explanation lines
#' @keywords internal
.format_resp_qualifier <- function(resps, multi_specs, n_total_resp) {
    resps <- resps[!is.na(resps) & resps != "" & resps != "All"]
    if (length(resps) == 0 || (!multi_specs && n_total_resp > 1)) {
        return("")
    }
    if (length(resps) == 1) {
        sprintf(" (%s)", resps[1])
    } else if (length(resps) <= 3) {
        sprintf(" (%s)", paste(resps, collapse = ", "))
    } else {
        sprintf(" (%s, +%d more)", paste(resps[1:2], collapse = ", "), length(resps) - 2)
    }
}

#' Format human-readable explanations of ROPE intervals
#'
#' @param resolved_df Data frame containing Response, Family, r_min, r_max
#' @param data Optional model data frame to check if Gaussian traits are standardized
#' @param standardized Logical. If TRUE, formats explanations on the standardized SD / latent liability scale.
#' @return List with 'lines' (character vector) and 'note' (optional string)
#' @keywords internal
explain_rope_spec <- function(resolved_df, data = NULL, standardized = FALSE) {
    if (is.null(resolved_df) || nrow(resolved_df) == 0) return(list(lines = character(0), note = NULL))

    # Keep all rows with non-NA values
    resolved_df <- unique(resolved_df[, c("Response", "Family", "r_min", "r_max")])

    total_responses <- unique(resolved_df$Response)
    n_total_resp <- length(total_responses)
    lines <- character(0)

    if (isTRUE(standardized)) {
        u_std <- unique(resolved_df[, c("r_min", "r_max"), drop = FALSE])
        for (i in seq_len(nrow(u_std))) {
            rmin <- u_std$r_min[i]
            rmax <- u_std$r_max[i]
            matching_rows <- resolved_df[resolved_df$r_min == rmin & resolved_df$r_max == rmax, , drop = FALSE]
            resps <- unique(matching_rows$Response)
            has_nongaus <- any(!tolower(as.character(matching_rows$Family)) %in% c("gaussian", "normal"))
            resp_qualifier <- .format_resp_qualifier(resps, multi_specs = nrow(u_std) > 1, n_total_resp = n_total_resp)
            unit_desc <- if (has_nongaus) "SD / latent liability SD" else "SD"
            lines <- c(lines, sprintf("Standardized%s: [%.2f, %.2f] (\u00b1%.2f %s)",
                                      resp_qualifier, rmin, rmax, max(abs(c(rmin, rmax))), unit_desc))
        }
        return(list(lines = lines, note = NULL))
    }

    # Find unique specs by Family, r_min, r_max
    u_specs <- unique(resolved_df[, c("Family", "r_min", "r_max")])
    unstandardized_vars <- character(0)

    for (i in seq_len(nrow(u_specs))) {
        fam  <- tolower(as.character(u_specs$Family[i]))
        rmin <- u_specs$r_min[i]
        rmax <- u_specs$r_max[i]

        matching_rows <- resolved_df[tolower(as.character(resolved_df$Family)) == fam &
                                     resolved_df$r_min == rmin &
                                     resolved_df$r_max == rmax, ]
        resps <- unique(matching_rows$Response)

        fam_title <- if (fam %in% c("gaussian", "normal")) "Gaussian"
                     else if (fam %in% c("poisson", "negbinomial", "zip", "zinb")) {
                         if (fam == "poisson") "Poisson" else toupper(fam)
                     } else {
                         paste0(toupper(substr(fam, 1, 1)), substr(fam, 2, nchar(fam)))
                     }

        resp_qualifier <- .format_resp_qualifier(resps, multi_specs = nrow(u_specs) > 1, n_total_resp = n_total_resp)

        if (fam %in% c("gaussian", "normal")) {
            lines <- c(lines, sprintf("%s%s: [%.2f, %.2f] (\u00b1%.2f SD if standardized)",
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
            lines <- c(lines, sprintf("%s%s: [%.2f, %.2f] (Odds Ratio [%.2f, %.2f], \u00b1%.2f latent SD)",
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
