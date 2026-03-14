#' Define JAGS Structure Implementation
#'
#' Modules implement this method to inject their specific JAGS code for covariance structures.
#' @param structure The structure object (e.g., phylo, list, custom class).
#' @param variable_name Name of the error variable (default "err").
#' @param ... Additional arguments.
#' @return A list containing `setup_code` (e.g. priors, matrix inversion) and `error_prior` (the likelihood/prior definition for residuals).
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
#' @export
nimble_family_optimization <- function(family, model_string, ...) {
    UseMethod("nimble_family_optimization")
}

#' Method for Structure Data Preparation (Matrix)
#' @export
prepare_structure_data.matrix <- function(
    structure,
    data,
    optimize = TRUE,
    ...
) {
    # Treat matrix as covariance, invert to get Precision
    # Precision matrix is expected by dmnorm
    return(list(data_list = list(Prec = solve(structure))))
}

#' Method for JAGS Structure Definition (Matrix)
#' @export
jags_structure_definition.matrix <- function(
    structure,
    variable_name = "err",
    optimize = TRUE,
    ...
) {
    # Matrix structures use dmnorm with the provided Precision matrix
    return(list(
        setup_code = NULL,
        error_prior = "dmnorm"
    ))
}

#' Default Method for D-Sep Transformation
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
#' @export
prepare_structure_data.default <- function(structure, data, ...) {
    return(list())
}

#' Default Method for Family Definition
#' @export
jags_family_definition.default <- function(family, response, predictors, ...) {
    return(NULL) # Triggers default Gaussian handling
}

#' Default Method for D-Sep Transformation
#' @export
transform_graph_for_dsep.default <- function(family, equations, ...) {
    return(equations) # No transformation
}

#' Default Method for NIMBLE Optimization
#' @export
nimble_family_optimization.default <- function(family, model_string, ...) {
    return(list(model_string = model_string, nimble_functions = list()))
}
