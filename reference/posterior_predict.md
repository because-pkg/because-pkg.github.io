# Posterior Predictive Samples

Generic function for generating posterior predictive draws from a fitted
model. For `because` model objects, see `posterior_predict.because`,
which returns an `[ndraws x N_obs]` matrix of simulated response values.
These draws are the basis for
[`pp_check`](https://because-pkg.github.io/because/reference/pp_check.md).

This function generates simulated data (\\y\_{rep}\\) from the posterior
distribution of a fitted `because` model. These draws are used for
posterior predictive checks (PPC) and model validation.

## Usage

``` r
posterior_predict(object, ...)

# S3 method for class 'because'
posterior_predict(
  object,
  resp = NULL,
  newdata = NULL,
  ndraws = NULL,
  re_formula = NULL,
  ...
)
```

## Arguments

- object:

  A `because` fit object.

- ...:

  Additional arguments (currently ignored).

- resp:

  Character string; the name of the response variable to predict. If
  `NULL`, takes the first response variable in the model.

- newdata:

  Optional new data frame for prediction. Currently ignored.

- ndraws:

  Integer; number of posterior draws to use. Defaults to all draws.

- re_formula:

  Formula or `NA`; determines which random effects to include. If `NULL`
  (default), all random effects are included (conditional prediction).
  If `NA`, no random effects are included (marginal prediction).

## Value

A matrix of posterior predictive draws.

A matrix of dimensions `[ndraws x N_obs]` containing simulated response
values.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(Y = rnorm(100), X = rnorm(100))
fit <- because(list(Y ~ X), data = df)

# Draw 200 posterior predictive samples for Y
yrep <- posterior_predict(fit, resp = "Y", ndraws = 200)
dim(yrep)  # [200 x 100]

# Marginal prediction (no random effects)
yrep_marg <- posterior_predict(fit, resp = "Y", ndraws = 200, re_formula = NA)
} # }
```
