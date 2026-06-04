# Causal Interventions for Because Models

This method implements Pearl's do-operator for Bayesian Structural
Equation Models fitted with `because`. It topologically propagates
interventions downstream through the causal graph to simulate
counterfactual outcomes.

## Usage

``` r
# S3 method for class 'because'
do(object, ..., ndraws = NULL, re_formula = NULL, raw_scale = FALSE)
```

## Arguments

- object:

  A fitted `because` model object.

- ...:

  Interventions specified as name-value pairs. Supported types include:

  - **Atomic (Hard)**: `temp = 25` (fixes value exactly)

  - **Shift (Additive)**: `temp = ~ . + 1.5` (adds 1.5 to natural
    values)

  - **Percentage Shift**: `temp = "+5%"` or `temp = "-10%"` (multiplies
    by 1.05 or 0.90)

  - **Stochastic**: `temp = ~ rnorm(n, mean = . + 1.5, sd = 0.05)`

  Note: In formulas, `.` evaluates to the current/natural values, and
  `n` evaluates to the number of elements required for random generator
  functions.

- ndraws:

  Integer; number of posterior draws to simulate. Defaults to all draws.

- re_formula:

  Formula; determines which random effects to condition on. Defaults to
  `NULL` (conditional on estimated groups, e.g., historical sites). Set
  to `NA` for marginal predictions over new, unmeasured groups.

- raw_scale:

  Logical; if `TRUE`, the engine will automatically detect variables
  that were z-scored using
  [`scale()`](https://rdrr.io/r/base/scale.html) prior to model fitting.
  It will unscale the data, apply your intervention on the raw metric,
  and then re-scale it back to the z-metric before propagating the
  counterfactual. Defaults to `FALSE`. **Important Note**: Base R's
  `data.frame(X = scale(X))` will silently strip the scaling attributes
  required for this feature to work. To preserve them, you must assign
  the scaled column to an existing data frame, e.g., `df$X <- scale(X)`.

## Value

An object of class `because_counterfactual` (a list) containing
`[ndraws x N_obs]` matrices for each endogenous variable under the
intervention.

## Warning on Scaled Data

If you scaled your data (e.g., z-scoring) prior to fitting the model to
obtain comparable path coefficients, **interventions must be specified
on the scaled metric** unless you set `raw_scale = TRUE`. For example,
if Temperature was z-scored,
`do(fit, Temp = ~ . + 1.5, raw_scale = FALSE)` increases Temperature by
**1.5 standard deviations**, not 1.5 degrees.

Furthermore, applying **percentage shifts** (e.g., `Temp = "+5%"`) to
z-scored data will multiply the z-scores by 1.05. Because the mean of
z-scored data is 0, this stretches the variance but does *not* shift the
mean (since 0 \* 1.05 = 0). Always carefully consider the scale of your
variables when designing interventions, and utilize `raw_scale = TRUE`
if you wish to intervene on the original metric.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a simple mediation model
df <- data.frame(X = rnorm(100), M = rnorm(100), Y = rnorm(100))
fit <- because(list(M ~ X, Y ~ M + X), data = df)

# 1. Hard intervention: set X to exactly 15 everywhere
res_hard <- do(fit, X = 15)
mean(res_hard$Y)

# 2. Shift intervention: increase X by 5 from its baseline
res_shift <- do(fit, X = ~ . + 5)

# 3. Stochastic shift: increase X by 5, with variance
res_stoch <- do(fit, X = ~ rnorm(n, mean = . + 5, sd = 0.5))
} # }
```
