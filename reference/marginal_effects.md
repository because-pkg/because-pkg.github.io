# Calculate Marginal Effects for Because Fit

Estimates the Average Marginal Effect (AME) or Marginal Effect at the
Mean (MEM) for all paths in a model. This is especially useful for
comparing coefficients across different families (e.g., Gaussian vs.
Binomial) as it converts them to a common scale (expected change in the
response).

## Usage

``` r
marginal_effects(
  fit,
  at = NULL,
  prob = 0.95,
  samples = 1000,
  multinomial_probabilities = FALSE
)
```

## Arguments

- fit:

  A because fit object.

- at:

  Character or list. If NULL (default), calculates Average Marginal
  Effects (AME). If "mean", calculates Marginal Effects at the Mean
  (MEM).

- prob:

  Numeric; probability mass for the credible interval (default 0.95).

- samples:

  Integer. Number of posterior samples to use (default 1000 for
  precision and consistency across plots).

- multinomial_probabilities:

  Logical. If TRUE, returns granular probability shifts for each
  category of multinomial (unordered) responses instead of a single
  expected value shift. Default FALSE.

## Value

A data frame with marginal effects per path.
