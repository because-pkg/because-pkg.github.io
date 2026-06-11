if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "Estimate",
        "LowerCI",
        "Test",
        "UpperCI",
        "edge_id",
        "edge_label",
        "edge_type",
        "label_display",
        "model_label",
        "occ_species",
        "significant",
        "weight_abs",
        "xend",
        "xmax",
        "xmin",
        "y",
        "yend",
        "ymax",
        "ymin",
        "Estimate",
        "LowerCI",
        "Test",
        "UpperCI", # repeated for clarity if needed
        "name",
        "name",
        "type",
        "Label",
        "Lower",
        "Path",
        "Significant",
        "Upper",
        "curvature",
        "dx",
        "dy",
        "to",
        "returnType"
    ))
}

# Internal environment for storing dynamic functions (e.g. for NIMBLE)
# to avoid polluting the global environment and satisfy CRAN requirements.
.because_env <- new.env(parent = emptyenv())

.onAttach <- function(libname, pkgname) {
  # Toggle this to TRUE whenever you release a major update to the Python
  # backend and want to notify your R users to upgrade.
  notify_numpyro_update <- FALSE 
  
  if (notify_numpyro_update) {
    # Only notify users who actually have the reticulate environment set up
    if (requireNamespace("reticulate", quietly = TRUE)) {
      if (reticulate::virtualenv_exists("because_env") || reticulate::condaenv_exists("because_env")) {
        packageStartupMessage(
          "----------------------------------------------------------------------\n",
          " NOTE: A new version of the NumPyro backend (because_py) is available!\n",
          " Run `install_because_numpyro()` to upgrade and get the latest features.\n",
          "----------------------------------------------------------------------"
        )
      }
    }
  }
}
