#' Define JAGS Structure Implementation
#'
#' Modules implement this method to inject their specific JAGS code for covariance structures.
#' @param structure The structure object (e.g., phylo, list, custom class).
#' @param variable_name Name of the error variable (default "err").
#' @param ... Additional arguments.
#' @return A list containing `setup_code` (e.g. priors, matrix inversion) and `error_prior` (the likelihood/prior definition for residuals).
#' @keywords internal
#' @export
jags_structure_definition <- function(structure, variable_name = "err", ...) {
    UseMethod("jags_structure_definition")
}

#' Default Method for JAGS Structure Definition
#'
#' Implements a standard precision-matrix-based correlation (dmnorm).
#' @keywords internal
#' @export
jags_structure_definition.default <- function(
    structure,
    variable_name = "err",
    ...
) {
    args <- list(...)
    s_name <- args$s_name %||% "Struct"
    loop_bound <- args$loop_bound %||% "N"
    k_idx <- args$category_index # e.g. "k" or "2"
    is_multi <- args$is_multi %||% FALSE
    zeros_name <- args$zeros_name %||% "zeros"
    i_index <- args$i_index %||% "i"

    prec_name <- paste0("Prec_", s_name)
    err_var <- paste0("u_", variable_name, "_", s_name)
    raw_var <- paste0("u_std_", variable_name, "_", s_name)
    
    is_unified <- any(vapply(c("phylo", "spatial", "group"), function(u) grepl(tolower(u), tolower(s_name)), logical(1)))
    if (is_unified) {
        sig_var <- paste0("sigma_", s_name, "_", variable_name)
    } else {
        sig_var <- paste0("sigma_", variable_name, "_", s_name)
    }

    # Category suffix for parameter names if needed
    k_suffix <- if (!is.null(k_idx)) paste0("_", k_idx) else ""

    # Unique loop variable for this structure and response
    j_idx <- paste0("g_", variable_name, "_", s_name)

    # Indexing for the error variable definition
    raw_index <- if (!is.null(k_idx)) {
        paste0(raw_var, "[1:", loop_bound, ", ", k_idx, "]")
    } else {
        paste0(raw_var, "[1:", loop_bound, "]")
    }

    # Precision matrix indexing (multi-structure)
    prec_index <- if (is_multi) {
        paste0(
            prec_name,
            "[1:",
            loop_bound,
            ", 1:",
            loop_bound,
            ", ",
            (k_idx %||% "1"),
            "]"
        )
    } else {
        paste0(prec_name, "[1:", loop_bound, ", 1:", loop_bound, "]")
    }

    engine <- args$engine %||% "jags"
    tau_var <- if (is_unified) {
        paste0("tau_u_", s_name, "_", variable_name)
    } else {
        paste0("tau_u_", variable_name, "_", s_name)
    }

    # Standard priors & dmnorm definition
    model_lines <- c(
        paste0(
            "  ",
            raw_index,
            " ~ dmnorm(",
            zeros_name,
            "[1:",
            loop_bound,
            "], ",
            prec_index,
            ")"
        ),
        paste0("  for (", j_idx, " in 1:", loop_bound, ") {"),
        if (!is.null(k_idx)) {
            paste0(
                "    ",
                err_var,
                "[",
                j_idx,
                ", ",
                k_idx,
                "] <- ",
                raw_var,
                "[",
                j_idx,
                ", ",
                k_idx,
                "] * ",
                sig_var,
                k_suffix
            )
        } else {
            paste0(
                "    ",
                err_var,
                "[",
                j_idx,
                "] <- ",
                raw_var,
                "[",
                j_idx,
                "] * ",
                sig_var,
                k_suffix
            )
        },
        "  }"
    )

    # Term to add to linear predictor
    # This is where we use the hierarchical index bridge (e.g. site_idx_obs[i])
    term_str <- if (!is.null(k_idx)) {
        paste0(err_var, "[", i_index, ", ", k_idx, "]")
    } else {
        paste0(err_var, "[", i_index, "]")
    }

    if (inherits(structure, "phylo") || inherits(structure, "multiPhylo")) {
        stop("Object of class 'phylo' detected. Phylogenetic modeling requires the 'because.phybase' extension.\nPlease install it via: remotes::install_github('because-pkg/because.phybase')")
    }

    return(list(
        model_lines = model_lines,
        term = term_str
    ))
}

#' Prepare Structure Data
#'
#' Modules implement this to process the structure object into data for JAGS.
#' @param structure The structure object.
#' @param data The model data.
#' @param ... Additional arguments.
#' @return A named list of data to be passed to JAGS (e.g., list(VCV = ...)).
#' @keywords internal
#' @export
prepare_structure_data <- function(structure, data, ...) {
    UseMethod("prepare_structure_data")
}

#' Default Method for Prepare Structure Data
#' @keywords internal
#' @export
prepare_structure_data.default <- function(structure, data, ...) {
    if (inherits(structure, "phylo") || inherits(structure, "multiPhylo")) {
        stop("Object of class 'phylo' detected. Phylogenetic modeling requires the 'because.phybase' extension.\nPlease install it via: remotes::install_github('because-pkg/because.phybase')")
    }
    return(list())
}

#' Define JAGS Family Implementation
#'
#' @param family The S3 family object (e.g., because_family_occupancy).
#' @param response The name of the response variable.
#' @param predictors Character vector of predictor names.
#' @param ... Additional arguments.
#' @return A list containing `model_code` (the likelihood block) and `monitor_params`.
#' @keywords internal
#' @export
jags_family_definition <- function(family, response, predictors, ...) {
    UseMethod("jags_family_definition")
}

#' Default Method for JAGS Family Definition
#' @keywords internal
#' @export
jags_family_definition.default <- function(family, response, predictors, ...) {
    args <- list(...)
    prior <- args$prior
    
    # Define robust defaults (Weakly Informative)
    # SD=10 (Prec=0.01) is a safe "Universal" default that fixes convergence 
    # without biasing human-scale results.
    beta_prior  <- prior$beta  %||% "dnorm(0, 0.01)"
    alpha_prior <- prior$alpha %||% "dnorm(0, 0.01)"
    
    # ... rest of the generic code remains the same but uses these variables
    return(NULL) # This is a generic, specific families override below
}

#' Transform Graph for D-Separation
#'
#' Modules implement this to rename variables in the causal graph to match their latent structure.
#' @param family The S3 family object.
#' @param equations The list of structural formulas.
#' @param ... Additional arguments.
#' @return The transformed list of formulas.
#' @keywords internal
#' @export
transform_graph_for_dsep <- function(family, equations, ...) {
    UseMethod("transform_graph_for_dsep")
}

#' Default Method for Transform Graph
#' @keywords internal
#' @export
transform_graph_for_dsep.default <- function(family, equations, ...) {
    return(equations)
}

#' NIMBLE Family Optimization
#'
#' Modules implement this to provide specialized NIMBLE optimizations (e.g. marginalization).
#' @param family The S3 family object.
#' @param model_string The current model code string.
#' @param ... Additional arguments.
#' @return A list with \code{model_string} and \code{nimble_functions} (list of nimbleFunction).
#' @keywords internal
#' @export
nimble_family_optimization <- function(family, model_string, ...) {
    UseMethod("nimble_family_optimization")
}

#' Default Method for NIMBLE Optimization
#' @keywords internal
#' @export
nimble_family_optimization.default <- function(family, model_string, ...) {
    return(list(model_string = model_string, nimble_functions = list()))
}

#' Expanded D-Sep Equation Hook
#'
#' Modules implement this to add supporting equations (e.g. detection models)
#' to a d-separation test.
#' @param family The S3 family list.
#' @param equations The full list of model equations.
#' @param dsep_equations The current list of equations for the test.
#' @param ... Additional arguments.
#' @return A modified list of equations.
#' @keywords internal
#' @export
dsep_equations_hook <- function(family, equations, dsep_equations, ...) {
    UseMethod("dsep_equations_hook")
}

#' Default Method for D-Sep Equations Hook
#' @keywords internal
#' @export
dsep_equations_hook.default <- function(
    family,
    equations,
    dsep_equations,
    ...
) {
    return(dsep_equations)
}

#' D-Sep Tree Hook
#'
#' Modules implement this to decide if a tree should be passed to a specific d-separation test.
#' @param tree The structure object (e.g. phylo).
#' @param test_eq The current d-separation test equation.
#' @param hierarchical_info Hierarchical metadata.
#' @param levels Level mapping.
#' @param ... Additional arguments.
#' @return The tree object or NULL.
#' @keywords internal
#' @export
dsep_tree_hook <- function(tree, test_eq, hierarchical_info, levels, ...) {
    UseMethod("dsep_tree_hook")
}

#' Default Method for D-Sep Tree Hook
#'
#' D-separation independence tests are pure regression tests; they should not
#' inherit the full phylogenetic or spatial covariance structure from the main
#' model, which would cause JAGS compilation failures (dimension mismatches or
#' undefined loop indices). Returning NULL ensures each test is fitted as a
#' plain mixed-effects regression without structured covariance.
#' @keywords internal
#' @export
dsep_tree_hook.default <- function(
    tree,
    test_eq,
    hierarchical_info,
    levels,
    ...
) {
    # If no hierarchical info, we return NULL for safety (standard models)
    if (is.null(hierarchical_info)) return(NULL)
    
    # In hierarchical models, we MUST return the full structure list/object
    # so that the sub-model can correctly attach structured covariance 
    # (like phylogeny) to ALL relevant levels (Response or Random Groups).
    # e.g., Metabolic_Rate (Species response) AND Abundance (via Species random effect).
    return(tree)
}

#' D-Sep Potential Latent Hook
#'
#' Modules implement this to remove variables from the 'potential latents' list (e.g. occupancy states).
#' @param family The S3 family list.
#' @param potential_latents Character vector of variables suspected to be latent.
#' @param ... Additional arguments.
#' @return A modified character vector of potential latents.
#' @keywords internal
#' @export
dsep_potential_latent_hook <- function(family, potential_latents, ...) {
    UseMethod("dsep_potential_latent_hook")
}

#' Default Method for D-Sep Potential Latents
#' @keywords internal
#' @export
dsep_potential_latent_hook.default <- function(family, potential_latents, ...) {
    return(potential_latents)
}

#' D-Sep Test Translation Hook
#'
#' Modules implement this to translate test equations (e.g. psi_Species ~ X to Species ~ X).
#' @param family The S3 family list.
#' @param test_eq The d-separation test equation.
#' @param ... Additional arguments.
#' @return A modified test equation.
#' @keywords internal
#' @export
dsep_test_translation_hook <- function(family, test_eq, ...) {
    UseMethod("dsep_test_translation_hook")
}

#' Default Method for D-Sep Translation
#' @keywords internal
#' @export
dsep_test_translation_hook.default <- function(family, test_eq, ...) {
    return(test_eq)
}

#' Normalize Equations Hook
#'
#' Modules implement this to normalize specialized aliases in formulas (e.g. psi_Species -> Species).
#' @param family The S3 family list.
#' @param equations The list of structural formulas.
#' @param ... Additional arguments.
#' @return A modified list of formulas.
#' @keywords internal
#' @export
normalize_equations_hook <- function(family, equations, ...) {
    UseMethod("normalize_equations_hook")
}

#' Default Method for Normalize Equations Hook
#' @keywords internal
#' @export
normalize_equations_hook.default <- function(family, equations, ...) {
    return(equations)
}

#' Initialization Hook
#'
#' Modules implement this to provide specialized initial values (e.g. latent states for occupancy).
#' @param family The S3 family list.
#' @param data The model data.
#' @param ... Additional arguments.
#' @return A named list of initial values.
#' @keywords internal
#' @export
get_inits_hook <- function(family, data, ...) {
    UseMethod("get_inits_hook")
}

#' Default Method for Initialization Hook
#' @keywords internal
#' @export
get_inits_hook.default <- function(family, data, ...) {
    inits <- list()
    # If family is a named vector, iterate and provide defaults for count families
    if (is.character(family) && !is.null(names(family))) {
        for (v in names(family)) {
            fam_type <- family[[v]]
            if (fam_type %in% c("poisson", "negbinomial", "zip", "zinb")) {
                if (v %in% names(data)) {
                    y_vals <- as.numeric(data[[v]])
                    # Intercept (alpha) initialization
                    m_y <- mean(y_vals, na.rm = TRUE)
                    inits[[paste0("alpha_", v)]] <- log(max(0.1, m_y))
                }

                # Distribution specific parameters
                if (fam_type %in% c("negbinomial", "zinb")) {
                    inits[[paste0("r_", v)]] <- 1
                }
                if (fam_type %in% c("zip", "zinb")) {
                    inits[[paste0("psi_", v)]] <- 0.5
                }
            }
        }
    }
    return(inits)
}

#' Needs Zero Inflation Hook
#'
#' Modules implement this to signal if a variable needs a zero-inflated likelihood.
#' @param family The S3 family list.
#' @param variable_name Name of the response variable.
#' @param ... Additional arguments.
#' @return Logical.
#' @keywords internal
#' @export
needs_zero_inflation_hook <- function(family, variable_name, ...) {
    UseMethod("needs_zero_inflation_hook")
}

#' Default Method for Zero Inflation Hook
#' @keywords internal
#' @export
needs_zero_inflation_hook.default <- function(family, variable_name, ...) {
    return(FALSE)
}

#' Variability Type Hook
#'
#' Modules implement this to provide default variability metadata (e.g. 'reps' for occupancy).
#' @param family The S3 family list.
#' @param variable_name Name of the variable.
#' @param ... Additional arguments.
#' @return Character string ('se', 'reps', or NULL).
#' @keywords internal
#' @export
get_variability_type_hook <- function(family, variable_name, ...) {
    UseMethod("get_variability_type_hook")
}

#' Default Method for Variability Type Hook
#' @keywords internal
#' @export
get_variability_type_hook.default <- function(family, variable_name, ...) {
    return(NULL)
}

#' Tree Extraction Hook
#'
#' Modules implement this to extract the appropriate tree/structure (e.g. tip extraction from multiPhylo).
#' @param structure The structural object.
#' @param ... Additional arguments.
#' @return A structure object (e.g. phylo) or NULL.
#' @keywords internal
#' @export
get_tree_hook <- function(structure, ...) {
    UseMethod("get_tree_hook")
}

#' Default Method for Tree Extraction Hook
#' @keywords internal
#' @export
get_tree_hook.default <- function(structure, ...) {
    return(NULL)
}

#' Monitor Variables Hook
#'
#' Modules implement this to specify which parameters to monitor for a given response (e.g. z_Y for occupancy).
#' @param family The S3 family list.
#' @param variable_name Name of the response variable.
#' @param ... Additional arguments.
#' @return Character vector of parameter names to monitor.
#' @keywords internal
#' @export
get_monitor_vars_hook <- function(family, variable_name, ...) {
    UseMethod("get_monitor_vars_hook")
}

#' Default Method for Monitor Variables Hook
#' @keywords internal
#' @export
get_monitor_vars_hook.default <- function(family, variable_name, ...) {
    return(variable_name)
}

#' Order Labels Hook
#'
#' Modules implement this to extract labels (e.g. tip labels from phylo or rownames from matrix).
#' @param structure The structural object.
#' @param ... Additional arguments.
#' @return Character vector of labels or NULL.
#' @keywords internal
#' @export
get_order_labels_hook <- function(structure, ...) {
    UseMethod("get_order_labels_hook")
}

#' Default Method for Order Labels Hook
#' @keywords internal
#' @export
get_order_labels_hook.default <- function(structure, ...) {
    if (is.matrix(structure)) {
        return(rownames(structure))
    }
    return(NULL)
}

#' Latent Child Hook
#'
#' Modules implement this if a variable is modeled via a latent state (e.g. occupancy Y -> z_Y).
#' @param family The S3 family list.
#' @param variable_name Name of the variable.
#' @param ... Additional arguments.
#' @return Logical.
#' @keywords internal
#' @export
is_latent_child_hook <- function(family, variable_name, ...) {
    UseMethod("is_latent_child_hook")
}

#' Default Method for Latent Child Hook
#' @keywords internal
#' @export
is_latent_child_hook.default <- function(family, variable_name, ...) {
    return(FALSE)
}

#' Structure Name Hook
#'
#' Modules implement this to provide specialized names for structure precision matrices.
#' @param structure The structural object.
#' @param ... Additional arguments.
#' @return Character string.
#' @keywords internal
#' @export
get_structure_name_hook <- function(structure, ...) {
    UseMethod("get_structure_name_hook")
}

#' Default Method for Structure Name
#' @keywords internal
#' @export
get_structure_name_hook.default <- function(structure, ...) {
    return(NULL)
}

#' @title Plot Path Coefficients
#' @name plot_coef
#' @description
#' Generic for plotting path coefficients from a fitted model object.
#' For \code{because} model objects, see \code{\link{plot_coef.because}}, which
#' produces a caterpillar plot of raw or marginal-effects estimates with
#' customisable colour schemes and significance highlighting.
#' @param object A fitted model object.
#' @param ... Additional arguments passed to the method.
#' @export
plot_coef <- function(object, ...) {
    UseMethod("plot_coef")
}

#' Posterior Predictive Samples
#'
#' @description
#' Generic function for generating posterior predictive draws from a fitted model.
#' For \code{because} model objects, see \code{\link{posterior_predict.because}},
#' which returns an \code{[ndraws x N_obs]} matrix of simulated response values.
#' These draws are the basis for \code{\link{pp_check}}.
#' @param object A fitted model object.
#' @param ... Additional arguments passed to the method.
#' @return A matrix of posterior predictive draws.
#' @export
posterior_predict <- function(object, ...) {
    UseMethod("posterior_predict")
}

#' Posterior Predictive Checks
#'
#' @description
#' Generic function for posterior predictive checks of a fitted model.
#' For \code{because} model objects, see \code{\link{pp_check.because}},
#' which wraps \code{bayesplot} functions (density overlay, histogram, test
#' statistics) and supports conditional or marginal prediction via
#' \code{re_formula}.
#' @param object A fitted model object.
#' @param ... Additional arguments passed to the method.
#' @return A ggplot object.
#' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}

#' @title Plot D-Separation Tests
#' @name plot_dsep
#' @param object A fitted model object.
#' @param ... Additional arguments.
#' @export
plot_dsep <- function(object, ...) {
    UseMethod("plot_dsep")
}

#' @keywords internal
#' @export
prepare_structure_data.matrix <- function(structure, data, ...) {
    # Treat matrix as covariance, invert to get Precision
    return(list(data_list = list(Prec = solve(structure))))
}

#' @keywords internal
#' @export
jags_structure_definition.matrix <- function(
    structure,
    variable_name = "err",
    ...
) {
    # Matrix uses default precision-matrix logic
    jags_structure_definition.default(
        structure,
        variable_name = variable_name,
        ...
    )
}
