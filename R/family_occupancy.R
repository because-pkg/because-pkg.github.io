#' Get Family Object for S3 Dispatch
#'
#' Converts a family name string into a family class object for S3 dispatch.
#' This is the gatekeeper that ensures module packages are installed before
#' allowing specialized families to be used.
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

    # Check for occupancy - require module package
    if (family_name == "occupancy") {
        # Check for because.detection (actual name) or because.occupancy (alias/legacy)
        if (
            !requireNamespace("because.detection", quietly = TRUE) &&
                !requireNamespace("because.occupancy", quietly = TRUE)
        ) {
            stop(
                "Occupancy models require the 'because.detection' package.\n",
                "Install with:\n",
                "  remotes::install_github('because-pkg/because.occupancy')\n",
                "Or:\n",
                "  devtools::install('/path/to/because.occupancy')"
            )
        }
    }

    return(structure(
        list(name = family_name),
        class = c(paste0("because_family_", family_name), "because_family")
    ))
}

# NOTE: S3 methods for specific families (e.g., occupancy) are provided by
# their respective module packages (e.g., because.occupancy).
# The base package does NOT define these methods to avoid dispatch conflicts.
