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

A data frame with marginal effects per path. Each row contains:

- Response:

  Name of the response variable.

- Predictor:

  Name of the predictor variable.

- Category:

  For categorical predictors: the comparison category (vs. reference).
  `NA` for continuous predictors.

- Effect:

  Posterior mean of the marginal effect.

- Lower:

  Lower bound of the credible interval (at `(1-prob)/2`).

- Upper:

  Upper bound of the credible interval (at `1-(1-prob)/2`).

- Family:

  Distribution family of the response.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a mixed-family SEM
df <- data.frame(
  Y = rbinom(100, 1, 0.4),
  M = rnorm(100),
  X = rnorm(100)
)
fit <- because(list(M ~ X, Y ~ M + X), data = df, family = c(Y = "binomial"))

# Average Marginal Effects (default) — on the response scale for all families
me <- marginal_effects(fit)
print(me)

# Marginal Effects at the Mean
me_mem <- marginal_effects(fit, at = "mean")
print(me_mem)

# Narrow credible interval
me_90 <- marginal_effects(fit, prob = 0.90)
} # }
```
