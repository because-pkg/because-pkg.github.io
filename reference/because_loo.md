# Calculate LOO-CV for a Because Model

Calculates Leave-One-Out Cross-Validation using Pareto Smoothed
Importance Sampling (PSIS-LOO) for a fitted Because model.

## Usage

``` r
because_loo(model, ...)
```

## Arguments

- model:

  A fitted model object of class `"because"` returned by
  [`because`](https://because-pkg.github.io/because/reference/because.md).
  **Note**: If the model was not fitted with `WAIC = TRUE` (so `log_lik`
  is missing), this function will automatically refit the model (using a
  short MCMC run) to calculate the likelihoods.

- ...:

  Additional arguments passed to
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

## Value

A `loo` object containing:

- estimates:

  Table with ELPD (expected log pointwise predictive density), LOO-IC,
  and p_loo

- diagnostics:

  Pareto k diagnostic values for each observation

- pointwise:

  Pointwise contributions to LOO-IC

## Details

LOO-CV (Leave-One-Out Cross-Validation) uses Pareto Smoothed Importance
Sampling to approximate leave-one-out predictive performance without
refitting the model N times. This is particularly useful for:

- Model comparison when models have different numbers of latent
  variables

- Identifying influential observations (via Pareto k diagnostics)

- Robust predictive performance assessment

**Pareto k diagnostics**:

- k \< 0.5: Excellent (all estimates reliable)

- 0.5 \< k \< 0.7: Good (estimates okay)

- 0.7 \< k \< 1: Problematic (estimates unreliable)

- k \> 1: Very problematic (refit model excluding these observations)

**Note on implementation**: This function extracts the pointwise
log-likelihoods calculated during model sampling (monitored as
`log_lik[i]`) when `WAIC = TRUE`. It does not re-compute likelihoods
from posterior samples in R, ensuring consistency with the fitted model
structure.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- because(data, tree, equations, WAIC = TRUE)
  loo_result <- because_loo(fit)
  print(loo_result)

  # Check for problematic observations
  plot(loo_result)

  # Compare models
  loo_compare(loo_result1, loo_result2)
} # }
```
