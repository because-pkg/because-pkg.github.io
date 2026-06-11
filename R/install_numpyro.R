#' Install NumPyro and Python dependencies for because
#'
#' This helper function installs the necessary Python dependencies (`numpyro`, `jax`) 
#' and the companion Python package (`because_py`) required to run the `numpyro` engine in `because()`.
#'
#' @param envname The name, or full path, of the Python environment in which the Python packages 
#'   are to be installed. If \code{NULL} (the default), it automatically detects whether RStudio 
#'   or \code{reticulate} has forced a specific environment (e.g. via \code{VIRTUAL_ENV}) and 
#'   uses that. Otherwise, it defaults to `"r-reticulate"`.
#' @param method Installation method. By default, `"auto"` automatically finds a method that will work 
#'   in the local environment. Change the default to force a specific installation method 
#'   (e.g., `"virtualenv"`, `"conda"`, or `"pip"`).
#' @param because_py_url The GitHub repository URL for the `because_py` Python module. 
#'   Update this to the official repository link once hosted.
#' @param ... Additional arguments passed to \code{reticulate::py_install()}.
#'
#' @export
install_because_numpyro <- function(envname = NULL, 
                                    method = "auto",
                                    because_py_url = "https://github.com/because-pkg/because_py/archive/main.zip",
                                    ...) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required to install Python dependencies. Please install it with: install.packages('reticulate')")
  }
  
  if (is.null(envname)) {
    if (Sys.getenv("VIRTUAL_ENV") != "") {
      envname <- Sys.getenv("VIRTUAL_ENV")
    } else if (Sys.getenv("RETICULATE_PYTHON") != "") {
      envname <- dirname(dirname(Sys.getenv("RETICULATE_PYTHON")))
    } else if (reticulate::py_available(initialize = FALSE)) {
      envname <- dirname(dirname(reticulate::py_config()$python))
    } else {
      envname <- "because_env"
    }
  }
  
  message("Preparing to install Python dependencies for the 'because' numpyro engine...")
  message("Target environment: ", envname)
  
  # The required packages: jax, numpyro, and the because_py github repo
  pkgs <- c("numpyro", "jax", "jaxlib", "networkx", "funsor", because_py_url)
  
  tryCatch({
    reticulate::py_install(
      packages = pkgs,
      envname = envname,
      method = method,
      pip = TRUE,
      ...
    )
    
    message("\nInstallation successful!")
    message("------------------------")
    message("You can now use the NumPyro engine in your models:")
    message('fit <- because(equations, data, engine = "numpyro")')
    message("\nNote: If you encounter issues finding the environment in a fresh R session, you may need to explicitly load it before calling because():")
    message(sprintf('reticulate::use_virtualenv("%s", required = TRUE)', envname))
    message("If you are using Miniconda/Anaconda, use this instead:")
    message(sprintf('reticulate::use_condaenv("%s", required = TRUE)', envname))
    
  }, error = function(e) {
    message("\nInstallation failed. Error details:")
    message(e$message)
    message("\nIf you do not have Python installed, reticulate may prompt you to install Miniconda first.")
    message("You can install it manually by running: reticulate::install_miniconda()")
  })
}
