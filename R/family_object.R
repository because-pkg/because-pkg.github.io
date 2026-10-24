#' Get Family Object for S3 Dispatch
#'
#' Converts a family name string into a family class object for S3 dispatch.
#' This provides a general mechanism for family extensions to be used.
#'
#' @param family_name The name of the family (e.g., "gaussian", "binomial", or extension families).
#'
#' @return An object of class `because_family_<name>` and `because_family`.
#'
#' @keywords internal
get_family_object <- function(family_name) {
    if (is.null(family_name) || is.na(family_name)) {
        return(structure(
            list(),
            class = c("because_family_gaussian", "because_family")
        ))
    }

    return(structure(
        list(name = family_name),
        class = c(paste0("because_family_", family_name), "because_family")
    ))
}

# NOTE: S3 methods for specialized domain families (e.g., detection models,
# capture-recapture, etc.) are provided by their respective extension packages.
# The base package defines only standard core families and generic S3 hooks.
