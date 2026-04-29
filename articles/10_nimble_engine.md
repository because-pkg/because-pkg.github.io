# High-Performance SEM with the NIMBLE Engine

## Introduction

While **because** uses **JAGS** as its default inference engine, it also
provides a high-performance alternative through the **NIMBLE** package
(de Valpine et al., 2017). NIMBLE compiles your BUGS models into C++ and
allows for advanced optimizations like parallel execution and adaptive
slice sampling.

This vignette explains how to leverage the NIMBLE engine to speed up
complex Structural Equation Models, particularly those involving
high-dimensional random effects or large datasets.

## Installation

The **NIMBLE** engine is an optional dependency for the **because**
package. Since it compiles models into C++, it requires a C++ compiler
to be installed on your system (e.g., **Rtools** on Windows, **Xcode
Command Line Tools** on macOS, or **build-essential** on Linux).

You can install NIMBLE directly from CRAN:

``` r
install.packages("nimble")
```

For detailed installation instructions and system-specific requirements,
please visit the [NIMBLE Download Page](https://r-nimble.org/download).

## Getting Started with NIMBLE

To use the NIMBLE engine, you simply need to set the `engine` argument
in the
[`because()`](https://because-pkg.github.io/because/reference/because.md)
function:

``` r
library(because)

# Define your SEM
equations <- list(
    Mass ~ Age + Sex,
    Lifespan ~ Mass + Habitat
)

# Dummy data for demonstration
my_data <- data.frame(
  Mass = rnorm(100),
  Age = rnorm(100),
  Sex = rbinom(100, 1, 0.5),
  Lifespan = rnorm(100),
  Habitat = factor(sample(c("Forest", "Grassland"), 100, replace = TRUE))
)

# Fit using NIMBLE instead of JAGS
fit <- because(
    equations = equations,
    data = my_data,
    engine = "nimble"
)
```

### Why use NIMBLE?

1.  **Speed**: NIMBLE compiles the model into C++ code, which can be
    significantly faster for complex models or large datasets.
2.  **Parallelization**: You can run multiple MCMC chains in parallel
    across your CPU cores.
3.  **Adaptive Sampling**: NIMBLE provides specialized samplers like
    `AF_slice` (Automated Factor Slice Sampler) that often mix better
    for correlated parameters than JAGS’s defaults.

## Parallel Execution

One of the easiest ways to speed up your analysis in **because** is to
run MCMC chains in parallel. While this feature is available for both
JAGS and NIMBLE engines, it is particularly beneficial when using
NIMBLE.

Since NIMBLE must compile the C++ code for each chain, running them in
parallel allows this compilation to happen simultaneously on different
CPU cores, significantly reducing the total “wait time” before sampling
begins.

``` r
fit <- because(
    equations = equations,
    data = my_data,
    engine = "nimble", # Works for "jags" too!
    parallel = TRUE,
    n.chains = 4,
    n.cores = 4
)
```

## Advanced Sampler Customization

NIMBLE allows you to control which sampling algorithm is used for every
node in your model. While `because` sets high-performance defaults (like
`AF_slice` for random effects and `slice` for precision parameters), you
can manually override these using the `nimble_samplers` argument.

``` r
# Override the default RW-MH sampler for a specific coefficient
fit <- because(
    equations = equations,
    data = my_data,
    engine = "nimble",
    nimble_samplers = list(
        "beta_Lifespan_Mass" = "slice"
    )
)
```

### When to Override Samplers?

While the default samplers are chosen for robustness, no single
algorithm is perfect for every dataset. You should consider using
`nimble_samplers` if you observe:

1.  **Low Effective Sample Size (ESS)**: One or two parameters have very
    few independent samples compared to the rest of the model.
2.  **High R-hat ($> 1.05$)**: A parameter fails to converge even with a
    long burn-in and many iterations.
3.  **“Fuzzy” Traceplots**: The traceplot looks like a “snake” or has
    long-term trends rather than looking like a consistent
    “caterpillar.”

**Pro-Tip**: If a specific regression coefficient ($\beta$) is mixing
poorly, swapping the default `RW` for a `slice` sampler is often the
most effective fix.

### Common Sampler Types in NIMBLE

When using the `nimble_samplers` argument, you can choose from several
built-in sampling algorithms. The most common ones for `because` models
are:

- **`slice`**: A robust scalar slice sampler. It is excellent for
  parameters with complex or non-standard distribution shapes and often
  mixes better than Random Walk MH.
- **`AF_slice`**: The “Automated Factor Slice Sampler.” This is a
  multivariate slice sampler that automatically detects and handles
  correlations between parameters. It is highly recommended for
  structured random effect vectors.
- **`RW`**: A standard scalar Random Walk Metropolis-Hastings sampler.
  It is very fast per iteration but might require more tuning (or more
  iterations) to reach a high effective sample size compared to slice
  sampling.
- **`RW_block`**: A multivariate Random Walk sampler. It samples a block
  of parameters jointly, which can be faster than individual `RW`
  samplers if the parameters are strongly correlated.
- **`categorical`**: A specialized sampler for discrete choice models
  (e.g., Multinomial and Ordinal models). It is extremely efficient for
  these data types.

### Default Optimizations and Sampler Selection

By default, the NIMBLE engine in `because` applies several technical
optimizations to match or exceed JAGS’s mixing quality. The following
table documents how samplers are assigned when using
`engine = "nimble"`:

| Parameter Type                                | NIMBLE Sampler (default)        | JAGS Sampler      | Why?                                                                                        |
|:----------------------------------------------|:--------------------------------|:------------------|:--------------------------------------------------------------------------------------------|
| **Linear coefficients** (`alpha_*`, `beta_*`) | **Slice** (if non-Gaussian)     | **Slice**         | Swapped from default `RW` to **Slice** for Poisson/Binomial to handle steep link curvature. |
| **Structured REs** (`u_std_.*[1:N]`)          | **Adaptive Slice** (`AF_slice`) | **Slice**         | Swapped from default `RW_block` to match JAGS’s robust mixing for correlated vectors.       |
| **Precision parameters** (`tau_u_.*`)         | **Slice** (`slice`)             | **Slice**         | Swapped from default `RW` to ensure stable convergence for variance components.             |
| **Dispersion / Inflation** (`r_*`, `psi_*`)   | **Slice** (`slice`)             | **Slice**         | Swapped from default `RW` to match JAGS’s robust mixing and handle boundary constraints.    |
| **Multinomial / Ordinal**                     | **Categorical** (`categorical`) | **Gibbs / Slice** | Specialized discrete sampler. Extremely fast and efficient.                                 |
| **Binomial / Poisson**                        | **Slice** (hardened)            | **Slice**         | Automatic hardening swaps the default `RW` for `slice` to prevent divergent transitions.    |

## Automatic Initialization and Stability

NIMBLE models can be sensitive to starting values, especially when using
non-linear link functions like
[`exp()`](https://rdrr.io/r/base/Log.html) (Poisson) or `logit()`
(Binomial). To ensure robust convergence, **because** implements several
automatic stabilization features:

### 1. Intelligent Intercept Starts

For count and binary responses, **because** automatically calculates the
mean of the observed data and initializes the corresponding intercept
(`alpha_Variable`) at the appropriate scale (e.g., `log(mean(y))`). This
ensures that chains begin in an identifiable region of the likelihood
rather than starting at `0.0`, which can lead to divergent transitions.

### 2. Stochastic Chain Jitter

To ensure that the `Rhat` statistic can accurately detect
non-convergence, **because** applies a stochastic jitter to the initial
values for every chain. Each chain starts at a slightly different point
in the parameter space, allowing the MCMC algorithm to explore the
posterior from multiple perspectives and escape locally flat regions.

### 3. Sampler Hardening

When **because** detects a non-Gaussian response, it automatically
“hardens” the sampler for the associated regression coefficients. While
NIMBLE defaults to a Random Walk (RW) sampler for most fixed effects,
`because` swaps this for a **Slice Sampler** for Poisson, Binomial, and
Negative Binomial models. This provides significantly more stability in
the presence of the extreme curvature characteristic of these link
functions.

### References

de Valpine, P., Turek, D., Paciorek, C. J., Anderson-Bergman, C., Lang,
D. T., & Bodik, R. (2017). Programming with models: Writing statistical
algorithms for general model structures with NIMBLE. *Journal of
Computational and Graphical Statistics*, 26(2), 403–413.
<https://doi.org/10.1080/10618600.2016.1172487>

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C.
(2021). Rank-normalization, folding, and localization: An improved R̂ for
assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667–718.
<https://doi.org/10.1214/20-BA1221>
