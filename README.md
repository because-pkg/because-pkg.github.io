# because

### Bayesian Estimation of Causal Effects

`because` provides a unified framework for specifying and fitting Bayesian structural equation models in `R` using [JAGS](http://mcmc-jags.sourceforge.net) or [NIMBLE](https://r-nimble.org).
The focus of `because` is on causal inference, providing tools to facilitate the correct estimation and testing of direct and indirect causal effects in complex systems and quantify uncertainty in the estimates.

> **Note on the logo**: The package hexagon sticker features the mathematical symbol for "because", represented as three dots in an inverted triangle ($\because$).

## Features

**because** simplifies the process of running complex Bayesian Structural Equation Models by automatically generating JAGS or NIMBLE code from standard R formulas. Key features include:

-   **Dual-Engine Support**: Choose between the simplicity of **JAGS** or the high-performance C++ backend of **NIMBLE**.
-   **Automatic Model Generation**: Builds complex BUGS/NIMBLE models directly from a list of structural equations.
-   **Generalized Covariance Structures**: Supports the specification of custom covariance structures in the data.
-   **Missing Data Support**: Imputes missing values in both response and predictor variables (assuming data is missing at random (MAR) or completely at random (MCAR).
-   **Measurement Error**: Accounts for measurement error providing repeated measures or known error variances.
-   **Distribution families**: Supports Gaussian, Binomial, Multinomial, Ordinal, Poisson,Negative Binomial as well as Zero Inflated Poisson (ZIP) and Zero Inflated Negative Binomial (ZINP) distributions.
-   **Categorical Predictors**: Automatic handling of factor variables with dummy variable expansion.
-   **Polynomial and Interaction Terms**: Easily include polynomial terms and interactions in structural equations.
-   **Hierachical Models**: Support for random effects and multi-level structures.
-   **Latent Variables**: Support for explicitly including Latent variables or modelling induced correlations from latent common causes with the Maximum Acyclic Graph (MAG) method.
-   **Causal inference tools**: Native support for d-separation and m-separation testing, mediation analysis, and **unified marginal effects** for comparing cross-family coefficients.
-   **Marginal Effects**: New `marginal_effects(fit)` function converts log-odds/latent scales to "Expected Change" units (0-1 probabilities or unit values) for easy comparison across different variable types.
-   **Parallel Computing**: Run MCMC chains in parallel on multi-core systems for faster computation.
-   **Visualisation tools**: Functions for visualizing model structures and posterior distributions.
-   **Extension packages**: Full S3 support for extension packages to build on the core functionality. 

## Installation

To install the **stable release** (`v1.2.2`), run:

``` r
remotes::install_github("because-pkg/because@v1.2.2", build_vignettes = TRUE)
```

To install the **latest development version** (unstable), run:

``` r
remotes::install_github("because-pkg/because", build_vignettes = TRUE)
```

## Prerequisites

Before using `because`, you need to have an inference engine installed.

### JAGS (Primary Engine)
-   **macOS**: `brew install jags` or download from [SourceForge](http://mcmc-jags.sourceforge.net).
-   **Windows**: Download installer from [SourceForge](http://mcmc-jags.sourceforge.net).
-   **Linux**: `sudo apt-get install jags`.

### NIMBLE (High-Performance Engine)
NIMBLE is an optional backend that requires a C++ compiler to be installed on your system:
-   **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
-   **macOS**: Install **Xcode Command Line Tools** via terminal: `xcode-select --install`.
-   **Linux**: Install `build-essential` (on Ubuntu/Debian) or equivalent development tools for your distribution.
-   **R Package**: `install.packages("nimble")`.
