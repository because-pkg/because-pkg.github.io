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

# Update these two variables whenever you release a new because_py backend!
.target_numpyro_version <- "0.1.1" 
.notify_numpyro_update <- FALSE

.onAttach <- function(libname, pkgname) {
  # Print the current R package version
  pkg_version <- utils::packageVersion(pkgname)
  packageStartupMessage(sprintf("This is because v%s", pkg_version))
  
  if (.notify_numpyro_update) {
    # Only notify users who actually have the reticulate environment set up
    if (requireNamespace("reticulate", quietly = TRUE)) {
      if (reticulate::virtualenv_exists("because_env") || reticulate::condaenv_exists("because_env")) {
        
        # Check if the user already installed this specific update
        config_dir <- tools::R_user_dir("because", which = "config")
        receipt_file <- file.path(config_dir, paste0("numpyro_receipt_", .target_numpyro_version, ".txt"))
        
        if (!file.exists(receipt_file)) {
          packageStartupMessage(
            "----------------------------------------------------------------------\n",
            sprintf(" NOTE: NumPyro backend (because_py) v%s is available!\n", .target_numpyro_version),
            " Run `install_because_numpyro()` to upgrade and get the latest features.\n",
            "----------------------------------------------------------------------"
          )
        }
      }
    }
  }
}

# Internal helper called by install_because_numpyro() to silence the message
.write_numpyro_update_receipt <- function() {
  config_dir <- tools::R_user_dir("because", which = "config")
  if (!dir.exists(config_dir)) {
    dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
  }
  receipt_file <- file.path(config_dir, paste0("numpyro_receipt_", .target_numpyro_version, ".txt"))
  file.create(receipt_file, showWarnings = FALSE)
}
