# Getting Started with because

## Getting started

**because** is an R package designed to easily perform causal inference
with Bayesian Structural Equation Models in JAGS. The package integrates
the methods proposed by von Hardenberg & Gonzalez-Voyer (2025) to fit
Phylogenetic Bayesian Structural Equation Models (PhyBaSE) and extends
them to other types of covariance structures (eg. spatial
autocorrelation, genetic relatedness etc.).

**because** main features:

- [Causal Inference with
  d-Separation](https://because-pkg.github.io/because/articles/02_dseparation.md):
  Testing conditional independencies implied by your causal model.
- [Custom Priors and Mechanistic
  Constraints](https://because-pkg.github.io/because/articles/03_custom_priors.md):
  Incorporating prior knowledge and mechanistic constraints.
- [Mediation
  Analysis](https://because-pkg.github.io/because/articles/04_mediation.md):
  Decomposing effects into direct and indirect components.
- [Deterministic
  Nodes](https://because-pkg.github.io/because/articles/05_deterministic_nodes.md):
  Interactions, thresholds and mathematical transformations.
- [Non-Gaussian Distribution
  Families](https://because-pkg.github.io/because/articles/08_non_gaussian_families.md):
  Modeling non-Gaussian data.
  - Gaussian (continuous data)
  - Binomial (binary/proportion data)
  - Multinomial (unordered categorical data)
  - Ordinal (ordered categorical data)
  - Poisson (count data)
  - Negative Binomial (overdispersed count data)
  - Zero-inflated Poisson (ZIP) and Negative Binomial (ZINB)
- [Bayesian Missing Data
  Imputation](https://because-pkg.github.io/because/articles/09_missing_data.md):
  Handling missing data through simultaneous Bayesian imputation within
  the SEM.
- [Model Diagnostics and
  Comparison](https://because-pkg.github.io/because/articles/10_model_diagnostics_comparison.md):
  Checking convergence, posterior predictive checks, comparing models
  with WAIC and LOO-CV.
- [Latent Variables and the MAG
  Approach](https://because-pkg.github.io/because/articles/11_latent_variables_mag.md):
  Causal inference in the presence of latent variables using the MAG
  approach by Shipley & Douda (2021).
- [Creating Custom Distribution
  Families](https://because-pkg.github.io/because/articles/12_custom_families.md):
  Extending `because` with user-defined JAGS distributions.
- *Advanced Model Specifications:* [random
  effects](https://because-pkg.github.io/because/articles/07_multilevel_models.md)
  (lme4-style `(1|group)` notation, mixed models, nested designs),
  polynomial terms, categorical predictors, measurement error.
- *[Multilevel & Hierarchical
  Data](https://because-pkg.github.io/because/articles/07_multilevel_models.md):*
  Multi-level data with variables at different hierarchical levels and
  correct degrees of freedom in d-separation tests.
- *Phylogenetic Path Analysis*: Using the Phylogenetic Bayesian
  Structural Equation Model approach (PhyBaSE, von Hardenberg &
  Gonzalez-Voyer, 2025) through the `because.phybase` extension.

### Quick Start

#### Installation

Before using **because**, you need to have **JAGS** (Just Another Gibbs
Sampler) installed on your machine.

- **macOS**: `brew install jags` or download from
  [SourceForge](http://mcmc-jags.sourceforge.net).
- **Windows**: Download installer from
  [SourceForge](http://mcmc-jags.sourceforge.net).
- **Linux**: `sudo apt-get install jags`.

After installing JAGS, you can install `because` from GitHub.

To install the **stable release** (`v1.2.7`):

``` r
remotes::install_github("because-pkg/because@v1.2.7", build_vignettes = TRUE)
```

To install the **latest development version**:

``` r
remotes::install_github("because-pkg/because", build_vignettes = TRUE)
```

Finally, load the package:

``` r
library(because)
```

#### Your First Model

The main function in because is
[`because()`](https://because-pkg.github.io/because/reference/because.md),
which compiles and fits your specified Structural Equation Model in
JAGS. The function is very rich with functionalities allowing to model
complex models with different error structures and hierarchically
structured data. However, here, to show the basic workflow, we will use
because() to fit a simple linear model involving only two variables.
Having only two variables this is not a typical SEM being equivalent to
a simple linear regression which can not be used to infer causality, but
it serves to illustrate the basic usage of the package.

Let’s start simulating two variables X and Y where Y is correlated to X:

``` r
# set seed for reproducibility
set.seed(67)

# Simulate predictor
n <- 100
X <- rnorm(n = n, mean = 50, sd = 10)

# Generate response with the chosen intercept (alpha) and slope (beta)
alpha <- 20
beta <- 0.5
Y <- alpha + beta * X + rnorm(n, mean = 0, sd = 10)

# Combine into data frame
sim.dat <- data.frame(X, Y)

# Fit linear model with lm() function for comparison
summary(lm(Y ~ X, data = sim.dat))
```

Now we can fit the same model using because. We need to specify the
structural equations (in this case only one) using R’s formula syntax
and then call
[`because()`](https://because-pkg.github.io/because/reference/because.md)
passing the equation and the data frame:

``` r
# Define the equations
equations <- list(Y ~ X)

# Fit the model
fit <- because(
  equations = equations,
  data = sim.dat
)

# View results
summary(fit)
```

To check for convergence of the MCMC chains, you can look at the Rhat
values in the summary output (should be \< 1.1). You can also plot the
trace plots of the MCMC samples:

``` r
# Plot trace plots
plot(fit$samples)
```

![](figures/traceplot_1.png)

`fit$samples` is an MCMC object, so you can use all coda functions to
analyze and plot the MCMC samples.

You can see how
[`because()`](https://because-pkg.github.io/because/reference/because.md)
translates your model into JAGS syntax calling `fit$model` or, before
fitting it, using
[`because_model()`](https://because-pkg.github.io/because/reference/because_model.md):

``` r
# Generate JAGS model code
jags_model_code <- because_model(
  equations = equations
)

cat(jags_model_code$model)
```

When specifing the equations you can include multiple predictors as well
as factors,interaction terms and polynomial terms following the
conventional R formula syntax.

### Next: Causal Inference with D-Separation

One of the main features of because is the ability to test your causal
model’s fit to the data using d-separation tests. D-separation tests
evaluate whether the conditional independencies implied by your causal
model hold in the data (Shipley, 2016). If they do not, this suggests
that your model may be misspecified and that you may need to add or
remove paths. More details on d-separation tests, how to interpret them
and a full tutorial can be found in the [Causal inference with
d-separation](https://because-pkg.github.io/because/articles/02_dseparation.md)
vignette.
