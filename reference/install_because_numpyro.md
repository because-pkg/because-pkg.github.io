# Install NumPyro and Python dependencies for because

This helper function installs the necessary Python dependencies
(`numpyro`, `jax`) and the companion Python package (`because_py`)
required to run the `numpyro` engine in
[`because()`](https://because-pkg.github.io/because/reference/because.md).

## Usage

``` r
install_because_numpyro(
  envname = NULL,
  method = "auto",
  because_py_url = "https://github.com/because-pkg/because_py/archive/main.zip",
  ...
)
```

## Arguments

- envname:

  The name, or full path, of the Python environment in which the Python
  packages are to be installed. If `NULL` (the default), it
  automatically detects whether RStudio or `reticulate` has forced a
  specific environment (e.g. via `VIRTUAL_ENV`) and uses that.
  Otherwise, it defaults to `"r-reticulate"`.

- method:

  Installation method. By default, `"auto"` automatically finds a method
  that will work in the local environment. Change the default to force a
  specific installation method (e.g., `"virtualenv"`, `"conda"`, or
  `"pip"`).

- because_py_url:

  The GitHub repository URL for the `because_py` Python module. Update
  this to the official repository link once hosted.

- ...:

  Additional arguments passed to
  [`reticulate::py_install()`](https://rstudio.github.io/reticulate/reference/py_install.html).
