# because

### Bayesian Estimation of Causal Effects

`because` provides a unified framework for specifying and fitting
Bayesian structural equation models in `R` using
[JAGS](https://mcmc-jags.sourceforge.io), [NIMBLE](https://r-nimble.org)
or [NumPyro](https://num.pyro.ai/en/stable/). The focus of `because` is
on causal inference, providing tools to facilitate the correct
estimation and testing of direct and indirect causal effects in complex
systems and quantify uncertainty in the estimates.

> **Note on the logo**: The package hexagon sticker features the
> mathematical symbol for “because”, represented as three dots in an
> inverted triangle (\because).

## Features

**because** simplifies the process of running complex Bayesian
Structural Equation Models by automatically generating JAGS, NIMBLE or
NumPyro code from standard R formulas. Key features include:

- **Multiple-Engine Support**: Choose between the simplicity of
  **JAGS**, the C++ speed of **NIMBLE**, or the HMC/NUTS performance of
  **NumPyro**.
- **Automatic Model Generation**: Builds complex JAGS/NIMBLE/NumPyro
  models directly from a list of structural equations.
- **Generalized Covariance Structures**: Supports the specification of
  custom covariance structures in the data.
- **Missing Data Support**: Imputes missing values in both response and
  predictor variables (assuming data is missing at random (MAR) or
  completely at random (MCAR)).
- **Measurement Error**: Accounts for measurement error providing
  repeated measures or known error variances.
- **Distribution families**: Supports Gaussian, Binomial, Multinomial,
  Ordinal, Poisson, Negative Binomial as well as Zero Inflated Poisson
  (ZIP) and Zero Inflated Negative Binomial (ZINB) distributions.
- **Categorical Predictors**: Automatic handling of factor variables
  with dummy variable expansion.
- **Polynomial and Interaction Terms**: Easily include polynomial terms
  and interactions in structural equations.
- **Hierarchical Models**: Support for random effects and multi-level
  structures.
- **Latent Variables**: Support for explicitly including Latent
  variables or modelling induced correlations from latent common causes
  with the Maximum Acyclic Graph (MAG) method.
- **Causal inference tools**: Native support for d-separation and
  m-separation testing, mediation analysis, and **unified marginal
  effects** for comparing cross-family coefficients.
- **Counterfactual Simulations**: A fully featured
  [`do()`](https://because-pkg.github.io/because/reference/do.md)
  operator to simulate structural policy interventions (atomic, shifts,
  percentage, stochastic) using Pearl’s do-calculus.
- **Marginal Effects**: `marginal_effects(fit)` converts log-odds/latent
  scales to “Expected Change” units (0-1 probabilities or unit values)
  for easy comparison across different variable types.
- **Parallel Computing**: Run MCMC chains in parallel on multi-core
  systems for faster computation.
- **Visualisation tools**: Functions for visualizing model structures
  and posterior distributions.
- **Extension packages**: Full S3 support for extension packages to
  build on the core functionality.

## Installation

To install the **stable release** (`v1.10.0`), run:

``` r

remotes::install_github("because-pkg/because@v1.10.0")
```

To install the **latest development version**, run:

``` r

remotes::install_github("because-pkg/because")
```

All tutorials and vignettes are available on the **[because
documentation website](https://because-pkg.github.io)**. If you also
want to build the vignettes locally for offline viewing
(`browseVignettes("because")`), run:

``` r

remotes::install_github("because-pkg/because", dependencies = TRUE, build_vignettes = TRUE)
```

> **Stay updated:** To receive email notifications when new versions or
> features are released, click **Watch → Custom → Releases** at the top
> right of the [GitHub
> repository](https://github.com/because-pkg/because), or check the
> [Changelog](https://because-pkg.github.io/because/news/).

## Prerequisites

Before using `because`, you need to have an inference engine installed.

### JAGS (Primary Engine)

- **macOS**: Download the **JAGS 5.0.2** universal installer
  (`JAGS-5.0.2-universal.pkg`) from [SourceForge (JAGS
  5.x)](https://sourceforge.net/projects/mcmc-jags/files/JAGS/5.x/macOS/)
  (or install via `brew install jags`).
- **Windows**: Download the **JAGS 5.0.2** installer
  (`Rtools45/JAGS-5.0.2.exe` for R 4.5+) from [SourceForge (JAGS
  5.x)](https://sourceforge.net/projects/mcmc-jags/files/JAGS/5.x/Windows/).
- **Linux**: Install **JAGS 5.x** from [SourceForge (JAGS 5.x
  Source)](https://sourceforge.net/projects/mcmc-jags/files/JAGS/5.x/Source/)
  if your distribution’s repository only provides JAGS 4.3.2
  (`sudo apt-get install jags`).

> **Important (JAGS 5.x vs. JAGS 4.x):** Current `rjags` (`5-5`) on CRAN
> requires **JAGS 5.x** (`JAGS >= 5.0.0`).
>
> **macOS users:** CRAN does not yet publish macOS binaries for
> `rjags 5-5` (`install.packages("rjags")` currently pulls `4-18`).
> After installing **JAGS 5.0.2**, install `rjags 5-5` *before*
> installing `because` (and restart R afterward) from the JAGS
> maintainer’s repository (Apple Silicon, R 4.6) or from CRAN source
> (requires `xcode-select --install`):
>
> ``` r
>
> install.packages("rjags", repos = c("https://ku-awdc.github.io/drat", "https://cloud.r-project.org"))
> # Or build from CRAN source:
> # install.packages("rjags", type = "source")
> ```
>
> **JAGS 4.3.2 systems (e.g., Ubuntu/Debian `apt`):** If your system
> only provides **JAGS 4.3.2** and you do not want to compile JAGS 5
> from source, install the matching archived `rjags 4-18` package
> *before* installing `because`:
>
> ``` r
>
> install.packages("https://cran.r-project.org/src/contrib/Archive/rjags/rjags_4-18.tar.gz", repos = NULL, type = "source")
> ```

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
all required Python packages (`numpyro`, `jax`, `jaxlib`, and the
`because_py` companion module):

``` r

install_because_numpyro()
```

If you encounter any issues finding your Python environment (or if you
are using Miniconda/Conda), we highly recommend referring to the
**[because NumPyro (because_py) Installation
Guide](https://because-pkg.github.io/because_py/articles/installation.html)**
for comprehensive troubleshooting and dedicated environment setups.

## Citation

To cite `because` in publications, please use:

``` bibtex
@Manual{,
  title = {because: Bayesian Estimation of Causal Effects},
  author = {Achaz von Hardenberg},
  year = {2026},
  note = {R package version 1.10.0},
  url = {https://because-pkg.github.io},
}
```

In text: \> von Hardenberg, A. (2026). *because: Bayesian Estimation of
Causal Effects*. R package version 1.10.0,
<https://github.com/achazhardenberg/because>.
