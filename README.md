# because

### Bayesian Estimation of Causal Effects

`because` provides a unified framework for specifying and fitting Bayesian structural equation models in `R` using [JAGS](http://mcmc-jags.sourceforge.net).
The focus of `because` is on causal inference, providing tools to facilitate the correct estimation and testng of direct and indirect causal effects in complex systems and quantify uncertainty in the estimates.

> **Note on the logo**: The package hexagon sticker features the mathematical symbol for "because", represented as three dots in an inverted triangle ($\because$).

## Features

**because** simplifies the process of running complex Bayesian Structural Equation Models by automatically generating JAGS code from standard R formulas. Key features include:

-   **Automatic JAGS Code Generation**: Builds models directly from a list of structural equations.
-   **Generalized Covariance Structures**: Supports the specification of custom covariance structures in the data.
-   **Missing Data Support**: Imputes missing values in both response and predictor variables (assuming data is missing at random (MAR) or completely at random (MCAR).
-   **Measurement Error**: Accounts for measurement error providing repeated measures or known error variances.
-   **Distribution families**: Supports Gaussian, Binomial, Multinomial, Ordinal, Poisson,Negative Binomial as well as Zero Inflated Poisson (ZIP) and Zero Inflated Negative Binomial (ZINP) distributions.
-   **Categorical Predictors**: Automatic handling of factor variables with dummy variable expansion.
-   **Polynomial and Interaction Terms**: Easily include polynomial terms and interactions in structural equations.
-   **Hierachical Models**: Support for random effects and multi-level structures.
-   **Latent Variables**: Support for explicitly including Latent variables or modelling induced correlations from latent common causes with the Maximum Acyclic Graph (MAG) method.
-   **Causal inference tools**: native support for d-separation and m-separation testing and mediation analysis.
-   **Parallel Computing**: Run MCMC chains in parallel on multi-core systems for faster computation.
-   **Visualisation tools**: Functions for visualizing model structures and posterior distributions.
-   **Extension packages**: Full S3 support for extension packages to build on the core functionality. 

## Installation

To install the **stable release** (`v1.0.0`), run:

``` r
remotes::install_github("because-pkg/because@v0.9.9", build_vignettes = TRUE)
```

To install the **latest development version** (unstable), run:

``` r
remotes::install_github("because-pkg/because", build_vignettes = TRUE)
```

## Prerequisites

Before using `because`, you need to have **JAGS** (Just Another Gibbs Sampler) installed on your machine.

-   **macOS**: `brew install jags` or download from [SourceForge](http://mcmc-jags.sourceforge.net).
-   **Windows**: Download installer from [SourceForge](http://mcmc-jags.sourceforge.net).
-   **Linux**: `sudo apt-get install jags`.
