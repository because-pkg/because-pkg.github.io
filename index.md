# because

### Bayesian Estimation of Causal Effects

`because` provides a unified framework for specifying and fitting
Bayesian structural equation models in `R` using
[JAGS](https://mcmc-jags.sourceforge.io) or
[NIMBLE](https://r-nimble.org). The focus of `because` is on causal
inference, providing tools to facilitate the correct estimation and
testing of direct and indirect causal effects in complex systems and
quantify uncertainty in the estimates.

> **Note on the logo**: The package hexagon sticker features the
> mathematical symbol for “because”, represented as three dots in an
> inverted triangle ($`\because`$).

## Features

**because** simplifies the process of running complex Bayesian
Structural Equation Models by automatically generating JAGS or NIMBLE
code from standard R formulas. Key features include:

- **Multiple-Engine Support**: Choose between the simplicity of **JAGS**
  or o the high-performance of **NIMBLE**
- **Automatic Model Generation**: Builds complex JAGS/NIMBLE/NUMPYRO
  models directly from a list of structural equations.
- **Generalized Covariance Structures**: Supports the specification of
  custom covariance structures in the data.
- **Missing Data Support**: Imputes missing values in both response and
  predictor variables (assuming data is missing at random (MAR) or
  completely at random (MCAR).
- **Measurement Error**: Accounts for measurement error providing
  repeated measures or known error variances.
- **Distribution families**: Supports Gaussian, Binomial, Multinomial,
  Ordinal, Poisson,Negative Binomial as well as Zero Inflated Poisson
  (ZIP) and Zero Inflated Negative Binomial (ZINB) distributions.
- **Categorical Predictors**: Automatic handling of factor variables
  with dummy variable expansion.
- **Polynomial and Interaction Terms**: Easily include polynomial terms
  and interactions in structural equations.
- **Hierachical Models**: Support for random effects and multi-level
  structures.
- **Latent Variables**: Support for explicitly including Latent
  variables or modelling induced correlations from latent common causes
  with the Maximum Acyclic Graph (MAG) method.
- **Causal inference tools**: Native support for d-separation and
  m-separation testing, mediation analysis, and **unified marginal
  effects** for comparing cross-family coefficients.
- **Counterfactual Simulations**: A fully featured operator to simulate
  structural policy interventions (atomic, shifts, percentage,
  stochastic) using Pearl’s do-calculus.
- **Marginal Effects**: New `marginal_effects(fit)` function converts
  log-odds/latent scales to “Expected Change” units (0-1 probabilities
  or unit values) for easy comparison across different variable types.
- **Parallel Computing**: Run MCMC chains in parallel on multi-core
  systems for faster computation.
- **Visualisation tools**: Functions for visualizing model structures
  and posterior distributions.
- **Extension packages**: Full S3 support for extension packages to
  build on the core functionality.

## Installation

To install the **stable release** (`v1.2.9`), run:

``` r

remotes::install_github("because-pkg/because@v1.2.9", build_vignettes = TRUE)
```

To install the **latest development version** (unstable), run:

``` r

remotes::install_github("because-pkg/because", build_vignettes = TRUE)
```

## Prerequisites

Before using `because`, you need to have an inference engine installed.

### JAGS (Primary Engine)

- **macOS**: `brew install jags` or download from
  [SourceForge](https://mcmc-jags.sourceforge.io).
- **Windows**: Download installer from
  [SourceForge](https://mcmc-jags.sourceforge.io).
- **Linux**: `sudo apt-get install jags`.

### NIMBLE (High-Performance Engine)

NIMBLE is an optional backend that requires a C++ compiler to be
installed on your system: - **Windows**: Install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). - **macOS**:
Install **Xcode Command Line Tools** via terminal:
`xcode-select --install`. - **Linux**: Install `build-essential` (on
Ubuntu/Debian) or equivalent development tools for your distribution. -
**R Package**: `install.packages("nimble")`.

### NumPyro (Python-based Engine)

NumPyro is an optional backend built on JAX, recommended for large
datasets or GPU-accelerated sampling. It runs from R via the
`reticulate` package. NumPyro requires Python 3.8 or later and the
`reticulate` R package. Start by installing `reticulate` if you don’t
have it yet:

``` r

install.packages("reticulate")
```

If you do not already have Python installed on your system, you can skip
installing it separately: `reticulate` can install Miniconda, which
bundles Python and conda together:

``` r

reticulate::install_miniconda()  # skip if Python is already installed
```

**One-step setup**: `because` provides a built-in helper that installs
all required Python packages (`numpyro`, `jax`, `jaxlib`, `networkx`,
`funsor`, and the `because_py` companion module):

``` r

install_because_numpyro()
```

The function automatically detects your active Python environment —
including RStudio project-level virtual environments (`.venv`) and any
existing `reticulate` configuration — falling back to `"r-reticulate"`
if none is found.

If in a fresh R session `reticulate` picks up a different Python
environment than the one where the packages were installed, `because`
will throw an error like:

    Error: Failed to import python module 'because.api'.
    Python is currently running from: /path/to/wrong/python
    Ensure because_py is installed in this exact environment.
    Try running: install_because_numpyro()

In this case, point `reticulate` to the correct environment before
calling
[`because()`](https://because-pkg.github.io/because/reference/because.md):

``` r

reticulate::use_virtualenv("r-reticulate", required = TRUE)
```
