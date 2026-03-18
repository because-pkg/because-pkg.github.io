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

#' Define JAGS Family Implementation
#'
#' Modules implement this method to inject specific JAGS model code for response families.
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

#' Default Method for D-Sep Tree Hook
#' @keywords internal
#' @export
dsep_tree_hook.default <- function(
    tree,
    test_eq,
    hierarchical_info,
    levels,
    ...
) {
    return(tree)
}

#' Default Method for D-Sep Potential Latents
#' @keywords internal
#' @export
dsep_potential_latent_hook.default <- function(family, potential_latents, ...) {
    return(potential_latents)
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

#' Default Method for Initialization Hook
#' @keywords internal
#' @export
get_inits_hook.default <- function(family, data, ...) {
    return(list())
}

#' Default Method for Zero Inflation Hook
#' @keywords internal
#' @export
needs_zero_inflation_hook.default <- function(family, variable_name, ...) {
    return(FALSE)
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

#' Default Method for Latent Child Hook
#' @keywords internal
#' @export
is_latent_child_hook.default <- function(family, variable_name, ...) {
    return(FALSE)
}

#' Default Method for Structure Name
#' @keywords internal
#' @export
get_structure_name_hook.default <- function(structure, ...) {
    return(NULL)
}

#' @rdname plot_dsep
#' @export
plot_dsep <- function(object, ...) {
    UseMethod("plot_dsep")
}

#' Method for Structure Data Preparation (Matrix)
#' @keywords internal
#' @export
prepare_structure_data.matrix <- function(
    structure,
    data,
    ...
) {
    # Treat matrix as covariance, invert to get Precision
    # Precision matrix is expected by dmnorm
    return(list(data_list = list(Prec = solve(structure))))
}

#' Method for JAGS Structure Definition (Matrix)
#' @keywords internal
#' @export
jags_structure_definition.matrix <- function(
    structure,
    variable_name = "err",
    ...
) {
    # Matrix structures use dmnorm with the provided Precision matrix
    return(list(
        setup_code = NULL,
        error_prior = "dmnorm"
    ))
}

#' Default Method for JAGS Structure Definition
#' @keywords internal
#' @export
jags_structure_definition.default <- function(
    structure,
    variable_name = "err",
    ...
) {
    # If structure is NULL or unrecognizable, return NULL
    # This triggers the default independent error model in because_model.R
    return(NULL)
}

#' Default Method for Structure Data Preparation
#' @keywords internal
#' @export
prepare_structure_data.default <- function(structure, data, ...) {
    return(list())
}

#' Default Method for Family Definition
#' @keywords internal
#' @export
jags_family_definition.default <- function(family, response, predictors, ...) {
    return(NULL) # Triggers default Gaussian handling
}

#' Default Method for D-Sep Transformation
#' @keywords internal
#' @export
transform_graph_for_dsep.default <- function(family, equations, ...) {
    return(equations) # No transformation
}

#' Default Method for NIMBLE Optimization
#' @keywords internal
#' @export
nimble_family_optimization.default <- function(family, model_string, ...) {
    return(list(model_string = model_string, nimble_functions = list()))
}
