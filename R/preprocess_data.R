#' Preprocess data and arguments for because()
#'
#' @param data The input data frame or list
#' @param equations The equations list
#' @param distribution The deprecated distribution arg
#' @param family The family arg
#' @param structure The structure arg
#' @return A list containing cleaned variables
#' @noRd
preprocess_because_args <- function(data, equations, distribution, family, structure) {
  # Input validation
  if (is.null(data)) {
    stop("Argument 'data' must be provided.")
  }
  if (is.null(equations)) {
    stop("Argument 'equations' must be provided.")
  }

  # --- Clean Data: Ensure matrix (e.g. from scale()) is a data.frame ---
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }

  # --- Clean Data: Coerce 1D matrices (e.g., from scale()) to vectors ---
  scale_info <- list()
  if (is.data.frame(data)) {
    for (nm in names(data)) {
      if (!is.null(attr(data[[nm]], "scaled:center"))) {
        scale_info[[nm]] <- list(
          center = attr(data[[nm]], "scaled:center"),
          scale = attr(data[[nm]], "scaled:scale")
        )
      }
      if (is.matrix(data[[nm]]) && ncol(data[[nm]]) == 1) {
        data[[nm]] <- as.numeric(data[[nm]])
      }
    }
  } else if (is.list(data)) {
    for (df_name in names(data)) {
      if (is.data.frame(data[[df_name]])) {
        for (nm in names(data[[df_name]])) {
          if (!is.null(attr(data[[df_name]][[nm]], "scaled:center"))) {
            scale_info[[nm]] <- list(
              center = attr(data[[df_name]][[nm]], "scaled:center"),
              scale = attr(data[[df_name]][[nm]], "scaled:scale")
            )
          }
          if (is.matrix(data[[df_name]][[nm]]) && ncol(data[[df_name]][[nm]]) == 1) {
            data[[df_name]][[nm]] <- as.numeric(data[[df_name]][[nm]])
          }
        }
      }
    }
  }

  # --- Backward Compatibility: distribution -> family ---
  if (!is.null(distribution)) {
    warning(
      "Argument 'distribution' is deprecated and will be removed in future versions. Please use 'family' instead."
    )
    if (is.null(family)) {
      family <- distribution
    }
  }

  # --- Family Object Normalization (Custom Families Support) ---
  family_objects <- list()
  if (!is.null(family)) {
    if (inherits(family, "because_family")) {
      family_objects[["_default"]] <- family
      family <- setNames(family$family, "_default")
    } else if (is.list(family) && !is.null(names(family))) {
      for (nm in names(family)) {
        if (inherits(family[[nm]], "because_family")) {
          family_objects[[nm]] <- family[[nm]]
          family[[nm]] <- family[[nm]]$family
        }
      }
      family <- unlist(family)
    }
    
    supported_families <- c(
      "gaussian", "normal", "poisson", "negbinomial", "zip", "zinb", 
      "binomial", "bernoulli", "multinomial", "ordinal"
    )
    
    for (node_name in names(family)) {
      specified_dist <- tolower(family[[node_name]])
      if (!(specified_dist %in% supported_families)) {
        stop(
          sprintf(
            "Error: The distribution family '%s' specified for node '%s' is not supported in 'because'. Supported families are: %s.", 
            specified_dist, node_name, paste(supported_families, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  }

  family_obj <- family
  if (!is.null(family) && ("occupancy" %in% family || "cmr" %in% family)) {
    class(family_obj) <- c("because_family_occupancy", class(family_obj))
  }
  structure_obj <- structure
  if (!is.null(structure)) {
    class(structure_obj) <- c("because_structure", class(structure_obj))
  }

  equations <- normalize_equations_hook(family_obj, equations, data = data)

  return(list(
    data = data,
    scale_info = scale_info,
    family = family,
    family_objects = family_objects,
    family_obj = family_obj,
    structure_obj = structure_obj,
    equations = equations
  ))
}
