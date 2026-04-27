# Generate Posterior Predictive Draws for Because Models

This function generates simulated data (\\y\_{rep}\\) from the posterior
distribution of a fitted `because` model. These draws are used for
posterior predictive checks (PPC) and model validation.

## Usage

``` r
# S3 method for class 'because'
posterior_predict(object, resp = NULL, ndraws = NULL, re_formula = NULL, ...)
```

## Arguments

- object:

  A `because` fit object.

- resp:

  Character string; the name of the response variable to predict. If
  `NULL`, takes the first response variable in the model.

- ndraws:

  Integer; number of posterior draws to use. Defaults to all draws.

- re_formula:

  Formula or `NA`; determines which random effects to include. If `NULL`
  (default), all random effects are included (conditional prediction).
  If `NA`, no random effects are included (marginal prediction).

- ...:

  Additional arguments (currently ignored).

## Value

A matrix of dimensions `[ndraws x N_obs]` containing simulated response
values.
