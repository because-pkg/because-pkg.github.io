


# ============================================================================
# Extension Generics for Family Hooks
# ============================================================================
#'
#' Extension packages (e.g., `because.phybase`) implement this generic to
#' tell the model builder which equations should skip standalone likelihood
#' and residual variance generation. Base `because` still emits `mu_X[i]` and
#' `alpha_X` priors (needed by the extension's JAGS code).
#'
#' @param family A family object (S3 class `because_family_*`)
#' @param response Character; the response variable name
#' @param all_responses Character vector of all response variable names in the model
#' @param ... Additional arguments
#' @return A list: `list(skip_likelihood = FALSE, skip_variance = FALSE)`
#' @export
is_auxiliary_equation <- function(family, response, all_responses, ...) {
    UseMethod("is_auxiliary_equation")
}

#' @export
is_auxiliary_equation.default <- function(family, response, all_responses, ...) {
    list(skip_likelihood = FALSE, skip_variance = FALSE)
}

#' Get the latent predictor variable name for a family
#'
#' Extension packages implement this to substitute latent states when their
#' variable is used as a predictor (e.g., occupancy uses `z_Species[i]`
#' instead of `Species[i]` in the linear predictor of dependent equations).
#'
#' @param family A family object (S3 class `because_family_*`)
#' @param predictor_name Character; the predictor variable name
#' @param predictor_index Character; the JAGS-indexed predictor string (e.g. `"Species[i]"`)
#' @param ... Additional arguments
#' @return Character; the JAGS-indexed string to use in the linear predictor
#' @export
get_latent_predictor_name <- function(family, predictor_name, predictor_index, ...) {
    UseMethod("get_latent_predictor_name")
}

#' @export
get_latent_predictor_name.default <- function(family, predictor_name, predictor_index, ...) {
    predictor_index  # Default: use the predictor as-is
}


#' Create a Custom Covariance Structure for use with because
#'
#' This helper function allows users to easily create custom covariance structures
#' by providing just a function that computes the precision (or covariance) matrix.
#' The function handles all the S3 method registration automatically.
#'
#' @param precision_fn A function that takes your structure data and returns
#'   a precision matrix (N x N). The first return value should be the precision matrix.
#'   Additional list elements can be returned and will be passed to JAGS as data.
#' @param name Optional character string; the name of your structure (e.g., "spatial_knn").
#'   Used to create the internal S3 class and JAGS precision matrix name (\code{Prec_<name>}).
#'   If omitted, an internal name is generated automatically.
#' @param description Optional description of the structure for documentation.
#'
#' @return A constructor function that creates structure objects of your custom class.
#'
#' @details
#' This function creates:
#' \itemize{
#'   \item A constructor function (returned) for creating structure objects
#'   \item S3 methods for \code{jags_structure_definition} (generates JAGS code)
#'   \item S3 methods for \code{prepare_structure_data} (prepares data for JAGS)
#' }
#'
#' The precision function should return either:
#' \itemize{
#'   \item A single matrix (the precision matrix)
#'   \item A list with \code{Prec} (the precision matrix) and optionally other data
#' }
#'
#' \strong{Engine compatibility:} Full BUGS code generation (including 3D precision
#' arrays for Bayesian model averaging) is supported for \code{engine = "jags"} and
#' \code{engine = "nimble"}. For \code{engine = "numpyro"}, standard 2D precision
#' matrices are passed via fallback, whereas 3D arrays require JAGS or NIMBLE.
#'
#' @examples
#' \dontrun{
#' # Create a spatial distance-decay structure
#' spatial_decay <- because_structure(
#'   precision_fn = function(coords, decay_rate = 0.1) {
#'     dist_mat <- as.matrix(dist(coords))
#'     W <- exp(-decay_rate * dist_mat)
#'     diag(W) <- 0
#'     D <- diag(rowSums(W))
#'     D - 0.99 * W  # Precision matrix
#'   }
#' )
#'
#' # Use it
#' my_struct <- spatial_decay(coords = my_coords, decay_rate = 0.2)
#' fit <- because(Y ~ X, data = data, structure = my_struct)
#' }
#'
#' @export
because_structure <- function(precision_fn = NULL, name = NULL, description = NULL) {
    # Support legacy positional call: because_structure("name", precision_fn)
    if (is.character(precision_fn) && is.function(name)) {
        tmp <- precision_fn
        precision_fn <- name
        name <- tmp
    }

    if (!is.function(precision_fn)) {
        stop("'precision_fn' must be a function")
    }

    if (is.null(name)) {
        idx <- if (exists(".custom_struct_count", envir = .because_env, inherits = FALSE)) {
            get(".custom_struct_count", envir = .because_env, inherits = FALSE) + 1L
        } else {
            1L
        }
        assign(".custom_struct_count", idx, envir = .because_env)
        name <- paste0("custom_structure_", idx)
    } else if (!is.character(name) || length(name) != 1) {
        stop("'name' must be a single character string or NULL")
    }

    # Clean name for R naming conventions
    clean_name <- gsub("[^a-zA-Z0-9_]", "_", name)
    class_name <- clean_name
    prec_data_name <- paste0("Prec_", clean_name)

    # Create the constructor function
    constructor <- function(...) {
        # Call user's precision function
        result <- precision_fn(...)

        # Normalize result
        if (is.matrix(result) || is.array(result)) {
            prec_matrix <- result
            extra_data <- list()
        } else if (is.list(result)) {
            if ("Prec" %in% names(result)) {
                prec_matrix <- result$Prec
                extra_data <- result[names(result) != "Prec"]
            } else {
                prec_matrix <- result[[1]]
                extra_data <- if (length(result) > 1) result[-1] else list()
            }
        } else {
            stop(
                "precision_fn must return a matrix, array, or a list containing one"
            )
        }

        # Validate precision matrix (supporting 2D matrix or 3D array for BMA)
        dims <- dim(prec_matrix)
        if (is.null(dims) || (length(dims) != 2 && length(dims) != 3)) {
            stop("Precision must be a square matrix or a 3D array of square matrices")
        }

        # Check squareness of the spatial dimensions
        if (length(dims) == 2) {
            if (dims[1] != dims[2]) stop("Precision matrix must be square")
            n <- dims[1]
        } else {
            if (dims[2] != dims[3]) stop("Slices in 3D precision array must be square matrices")
            n <- dims[2]
        }

        # Create structure object
        structure(
            c(
                list(
                    Prec = prec_matrix,
                    n = n,
                    name = clean_name,
                    description = description
                ),
                extra_data
            ),
            class = c(
                class_name,
                "because_custom_structure",
                "because_structure"
            )
        )
    }

    # Register S3 methods in the calling environment
    # jags_structure_definition method
    jags_method <- function(
        structure,
        variable_name = paste0("u_", class_name),
        ...
    ) {
        setup_code <- c(
            paste0("    # Custom Structure: ", clean_name),
            if (!is.null(description)) paste0("    # ", description) else NULL
        )

        # Detect BMA (3D precision array)
        is_bma <- !is.null(dim(structure$Prec)) && length(dim(structure$Prec)) == 3

        if (is_bma) {
            n_slices <- dim(structure$Prec)[1]
            # Use variable_name (e.g. log_mass) instead of clean_name (ou_kernel)
            # to allow trait-specific K indices in the same model
            k_name <- paste0("K_", variable_name)
            p_name <- paste0("p_", variable_name)

            setup_code <- c(
                setup_code,
                paste0("    ", k_name, " ~ dcat(", p_name, "[1:", n_slices, "])"),
                paste0("    for (k in 1:", n_slices, ") { ", p_name, "[k] <- 1/", n_slices, " }")
            )

            # 3D Indexing for BMA using the unique K name
            error_prior <- paste0(
                "    ",
                variable_name,
                "[1:N] ~ dmnorm(zeros[1:N], tau_",
                clean_name,
                " * ",
                prec_data_name,
                "[", k_name, ", 1:N, 1:N])"
            )
        } else {
            # Standard multivariate normal with 2D precision matrix
            error_prior <- paste0(
                "    ",
                variable_name,
                "[1:N] ~ dmnorm(zeros[1:N], tau_",
                clean_name,
                " * ",
                prec_data_name,
                "[1:N, 1:N])"
            )
        }

        list(
            setup_code = setup_code,
            error_prior = error_prior
        )
    }

    # prepare_structure_data method
    prepare_method <- function(structure, data, optimize = TRUE, ...) {
        data_list <- list()
        data_list[[prec_data_name]] <- structure$Prec

        # Add any extra data from the structure
        for (nm in names(structure)) {
            if (!nm %in% c("Prec", "n", "name", "description", "class")) {
                data_list[[nm]] <- structure[[nm]]
            }
        }

        list(
            structure_object = structure,
            data_list = data_list
        )
    }

    # Assign methods to the global environment
    method_name_jags <- paste0("jags_structure_definition.", class_name)
    method_name_prep <- paste0("prepare_structure_data.", class_name)

    assign(method_name_jags, jags_method, envir = .because_env)
    assign(method_name_prep, prepare_method, envir = .because_env)

    # Register as S3 methods
    # We find the generic in the 'because' namespace (where they are defined)
    # This ensures that even in non-interactive sessions (Rscript),
    # the lookup doesn't fail.
    ns <- tryCatch(asNamespace("because"), error = function(e) .GlobalEnv)
    registerS3method(
        "jags_structure_definition",
        class_name,
        jags_method,
        envir = ns
    )
    registerS3method(
        "prepare_structure_data",
        class_name,
        prepare_method,
        envir = ns
    )

    message(
        "Created custom structure '",
        clean_name,
        "'\n",
        "  - Class: ",
        class_name,
        "\n",
        "  - Precision data: ",
        prec_data_name,
        "\n",
        "  - Registered S3 methods: jags_structure_definition.",
        class_name,
        ", prepare_structure_data.",
        class_name
    )

    return(constructor)
}


#' Create a Custom Family (Distribution) for use with because
#'
#' This helper function allows users to easily create custom distribution families
#' by providing the JAGS likelihood code. The function handles S3 method registration.
#'
#' @param jags_likelihood A character string with the JAGS likelihood code.
#'   Use placeholders: \code{\{response\}} for variable name, \code{\{mu\}} for mean,
#'   \code{\{i\}} for loop index, \code{\{tau\}} for precision parameter.
#' @param link Link function name (default: "identity"). Options: "identity", "log", "logit".
#' @param extra_priors Character vector of additional JAGS prior statements.
#'   Use \code{\{response\}} placeholder for variable-specific naming.
#' @param name Optional character string; the name of your family (e.g., "student_t").
#'   If omitted, an internal name is generated automatically.
#' @param description Optional description of the family.
#'
#' @return A constructor function for the family that can be passed to because(..., family = ...).
#'
#' @details
#' \strong{Engine compatibility:} Custom families defined via \code{because_family()}
#' inject BUGS-syntax likelihoods and priors into the model specification and are
#' supported for \code{engine = "jags"} and \code{engine = "nimble"} (provided the
#' BUGS distribution is available in NIMBLE). They are not supported for
#' \code{engine = "numpyro"}.
#'
#' @examples
#' \dontrun{
#' # Create a Student-t family for robust regression
#' student_t <- because_family(
#'   jags_likelihood = "{response}[{i}] ~ dt({mu}[{i}], {tau}, df_{response})",
#'   extra_priors = c("df_{response} ~ dunif(2, 100)")
#' )
#'
#' # Use it
#' fit <- because(Y ~ X, data = data, family = student_t())
#' }
#'
#' @export
because_family <- function(
    jags_likelihood = NULL,
    link = "identity",
    extra_priors = NULL,
    name = NULL,
    description = NULL
) {
    # Support legacy positional call: because_family("name", jags_likelihood, ...)
    if (
        is.null(name) &&
        is.character(jags_likelihood) && length(jags_likelihood) == 1 && !grepl("~", jags_likelihood) &&
        is.character(link) && length(link) >= 1 && any(grepl("~", link))
    ) {
        name <- jags_likelihood
        jags_likelihood <- link
        link <- "identity"
    }

    if (!is.character(jags_likelihood) || length(jags_likelihood) < 1) {
        stop("'jags_likelihood' must be a character string with JAGS code")
    }

    if (is.null(name)) {
        idx <- if (exists(".custom_family_count", envir = .because_env, inherits = FALSE)) {
            get(".custom_family_count", envir = .because_env, inherits = FALSE) + 1L
        } else {
            1L
        }
        assign(".custom_family_count", idx, envir = .because_env)
        dist_match <- regmatches(jags_likelihood[1], regexpr("~\\s*[a-zA-Z0-9_.]+", jags_likelihood[1]))
        dist_tok <- if (length(dist_match) == 1 && nchar(dist_match) > 1) {
            gsub("^~\\s*", "", dist_match)
        } else {
            "family"
        }
        name <- paste0("custom_", dist_tok, "_", idx)
    } else if (!is.character(name) || length(name) != 1) {
        stop("'name' must be a single character string or NULL")
    }

    clean_name <- gsub("[^a-zA-Z0-9_]", "_", name)
    class_name <- paste0("because_family_", clean_name)

    # Create the jags_family_likelihood method
    likelihood_method <- function(
        family,
        response,
        predictors = NULL,
        suffix = "",
        has_structure = FALSE,
        link_override = NULL,
        ...
    ) {
        # Use provided link or override
        use_link <- if (!is.null(link_override)) link_override else link

        # Process likelihood template
        likelihood <- jags_likelihood
        mu_var <- paste0("mu_", response, suffix)
        tau_var <- paste0("tau_res_", response, suffix)

        likelihood <- gsub("\\{response\\}", response, likelihood)
        likelihood <- gsub("\\{mu\\}", mu_var, likelihood)
        likelihood <- gsub("\\{tau\\}", tau_var, likelihood)
        likelihood <- gsub("\\{i\\}", "i", likelihood)

        # Add indentation
        likelihood_code <- paste0("    ", likelihood)

        # Process extra priors
        prior_code <- NULL
        if (!is.null(extra_priors)) {
            prior_code <- sapply(
                extra_priors,
                function(p) {
                    p <- gsub("\\{response\\}", response, p)
                    p <- gsub("\\{tau\\}", tau_var, p)
                    paste0("  ", p)
                },
                USE.NAMES = FALSE
            )
        }

        list(
            likelihood_code = likelihood_code,
            prior_code = prior_code,
            data_requirements = NULL
        )
    }

    # Register the S3 method
    method_name <- paste0("jags_family_likelihood.", class_name)
    assign(method_name, likelihood_method, envir = .because_env)

    ns <- tryCatch(asNamespace("because"), error = function(e) .GlobalEnv)
    registerS3method(
        "jags_family_likelihood",
        class_name,
        likelihood_method,
        envir = ns
    )

    # Create a constructor function that returns a family object
    constructor <- function() {
        structure(
            list(
                family = clean_name,
                link = link,
                description = description
            ),
            class = c(class_name, "because_custom_family", "because_family")
        )
    }

    message(
        "Created custom family '",
        clean_name,
        "'\n",
        "  - Class: ",
        class_name,
        "\n",
        "  - Link: ",
        link,
        "\n",
        "  - Registered S3 method: jags_family_likelihood.",
        class_name
    )

    return(constructor)
}
