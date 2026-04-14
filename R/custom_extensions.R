#' Create a Custom Covariance Structure for use with because
#'
#' This helper function allows users to easily create custom covariance structures
#' by providing just a function that computes the precision (or covariance) matrix.
#' The function handles all the S3 method registration automatically.
#'
#' @param name Character string; the name of your structure (e.g., "spatial_knn").
#'   This will be used to create the S3 class and methods.
#' @param precision_fn A function that takes your structure data and returns
#'   a precision matrix (N x N). The first return value should be the precision matrix.
#'   Additional list elements can be returned and will be passed to JAGS as data.
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
#' @examples
#' \dontrun{
#' # Create a spatial distance-decay structure
#' spatial_decay <- because_structure(
#'   name = "spatial_decay",
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
because_structure <- function(name, precision_fn, description = NULL) {
    # Validate inputs
    if (!is.character(name) || length(name) != 1) {
        stop("'name' must be a single character string")
    }
    if (!is.function(precision_fn)) {
        stop("'precision_fn' must be a function")
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
    ns <- asNamespace("because")
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
#' @param name Character string; the name of your family (e.g., "student_t").
#' @param jags_likelihood A character string with the JAGS likelihood code.
#'   Use placeholders: \code{\{response\}} for variable name, \code{\{mu\}} for mean,
#'   \code{\{i\}} for loop index, \code{\{tau\}} for precision parameter.
#' @param link Link function name (default: "identity"). Options: "identity", "log", "logit".
#' @param extra_priors Character vector of additional JAGS prior statements.
#'   Use \code{\{response\}} placeholder for variable-specific naming.
#' @param description Optional description of the family.
#'
#' @return A constructor function for the family that can be passed to because(..., family = ...).
#'
#' @examples
#' \dontrun{
#' # Create a Student-t family for robust regression
#' student_t <- because_family(
#'   name = "student_t",
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
    name,
    jags_likelihood,
    link = "identity",
    extra_priors = NULL,
    description = NULL
) {
    if (!is.character(name) || length(name) != 1) {
        stop("'name' must be a single character string")
    }
    if (!is.character(jags_likelihood) || length(jags_likelihood) < 1) {
        stop("'jags_likelihood' must be a character string with JAGS code")
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
        tau_var <- paste0("tau_e_", response, suffix)

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

    registerS3method(
        "jags_family_likelihood",
        class_name,
        likelihood_method
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
        class_name,
        "\n",
        "  - Usage: because(..., family = ",
        clean_name,
        "())"
    )

    return(constructor)
}
