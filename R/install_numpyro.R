#' Install NumPyro and Python dependencies for because
#'
#' This helper function installs the necessary Python dependencies (`numpyro`, `jax`) 
#' and the companion Python package (`because_py`) required to run the `numpyro` engine in `because()`.
#'
#' @param envname The name, or full path, of the Python environment in which the Python packages 
#'   are to be installed. Defaults to `"r-reticulate"`.
#' @param method Installation method. By default, `"auto"` automatically finds a method that will work 
#'   in the local environment. Change the default to force a specific installation method 
#'   (e.g., `"virtualenv"`, `"conda"`, or `"pip"`).
#' @param because_py_url The GitHub repository URL for the `because_py` Python module. 
#'   Update this to the official repository link once hosted.
#' @param ... Additional arguments passed to \code{reticulate::py_install()}.
#'
#' @export
install_because_numpyro <- function(envname = "r-reticulate", 
                                    method = "auto",
                                    because_py_url = "git+https://github.com/because-pkg/because_py.git",
                                    ...) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("The 'reticulate' package is required to install Python dependencies. Please install it with: install.packages('reticulate')")
  }
  
  message("Preparing to install Python dependencies for the 'because' numpyro engine...")
  message("Target environment: ", envname)
  
  # The required packages: jax, numpyro, and the because_py github repo
  pkgs <- c("jax", "jaxlib", "numpyro", because_py_url)
  
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
    
  }, error = function(e) {
    message("\nInstallation failed. Error details:")
    message(e$message)
    message("\nIf you do not have Python installed, reticulate may prompt you to install Miniconda first.")
    message("You can install it manually by running: reticulate::install_miniconda()")
  })
}
