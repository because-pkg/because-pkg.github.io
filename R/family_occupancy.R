#' Get Family Object for S3 Dispatch
#'
#' Converts a family name string into a family class object for S3 dispatch.
#' This provides a general mechanism for family extensions to be used.
#'
#' @param family_name The name of the family (e.g., "occupancy", "gaussian").
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

# NOTE: S3 methods for specific families (e.g., occupancy) are provided by
# their respective extension module packages.
# The base package does NOT define these methods to avoid dispatch conflicts.
