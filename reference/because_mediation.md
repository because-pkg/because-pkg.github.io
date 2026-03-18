# Automated Mediation Analysis

Decomposes the total effect of an exposure on an outcome into direct and
indirect pathways by tracing paths in the fitted DAG and multiplying
posterior coefficients.

## Usage

``` r
because_mediation(fit, exposure, outcome, prob = 0.95)
```

## Arguments

- fit:

  A fitted object from
  [`because()`](https://because-pkg.github.io/because/reference/because.md).

- exposure:

  Character string; the name of the exposure variable.

- outcome:

  Character string; the name of the outcome variable.

- prob:

  Numeric; probability mass for the credible interval (default 0.95).

## Value

A list containing:

- `paths`: A data frame summarizing each path (Path string, Mean, SD,
  CI).

- `summary`: A data frame summarizing Total, Direct, and Total Indirect
  effects.

- `samples`: A matrix of posterior samples for each path and the total
  effect.

## Details

This function reconstructs the causal graph from the `parameter_map`
stored in the fit object. It uses `igraph` to find all simple paths from
`exposure` to `outcome`. For each path, it identifies the corresponding
regression coefficients (`beta_...`) and computes their product across
all MCMC samples.

**Note:** This assumes linear relationships for indirect effects
(\\\beta_1 \times \beta_2\\). For non-linear models, this is an
approximation of the average causal effect.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming a dataset 'df' exists with variables X, M, and Y
equations <- list(M ~ X, Y ~ M + X)
fit <- because(equations, data = df)
med_results <- because_mediation(fit, exposure = "X", outcome = "Y")
print(med_results$summary)
} # }
```
