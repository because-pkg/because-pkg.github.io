# High-Performance SEM with the NIMBLE Engine

## Introduction

While **because** uses **JAGS** as its default inference engine, it also
provides a high-performance alternative through the **NIMBLE** package.
NIMBLE compiles your BUGS models into C++ and allows for advanced
optimizations like parallel execution and adaptive slice sampling.

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

| Parameter Type                                | NIMBLE Sampler (default)        | JAGS Sampler      | Why?                                                                                     |
|:----------------------------------------------|:--------------------------------|:------------------|:-----------------------------------------------------------------------------------------|
| **Linear coefficients** (`alpha_*`, `beta_*`) | **Gibbs (conjugate)**           | **Gibbs**         | Identical. Both use closed-form conjugate updates. Fast and exact.                       |
| **Residual precision** (`tau_e_*`)            | **Gibbs (conjugate)**           | **Gibbs**         | Gamma prior + Normal likelihood = conjugate.                                             |
| **Structured REs** (`u_std_.*[1:N]`)          | **Adaptive Slice** (`AF_slice`) | **Slice**         | Swapped from default `RW_block` to match JAGS’s robust mixing for correlated vectors.    |
| **Precision parameters** (`tau_u_.*`)         | **Slice** (`slice`)             | **Slice**         | Swapped from default `RW` to ensure stable convergence for variance components.          |
| **Dispersion / Inflation** (`r_*`, `psi_*`)   | **Slice** (`slice`)             | **Slice**         | Swapped from default `RW` to match JAGS’s robust mixing and handle boundary constraints. |
| **Multinomial / Ordinal**                     | **Categorical** (`categorical`) | **Gibbs / Slice** | Specialized discrete sampler. Extremely fast and efficient.                              |
| **Binomial / Poisson**                        | **RW-MH scalar** (`RW`)         | **Slice**         | `RW` is faster per iteration; parallelization compensates for mixing differences.        |
