# Posterior Predictive Checks for Because Fits

A wrapper around
[`bayesplot::ppc_dens_overlay`](https://mc-stan.org/bayesplot/reference/PPC-distributions.html)
and other PPC functions for `because` model objects.

## Usage

``` r
# S3 method for class 'because'
pp_check(
  object,
  resp = NULL,
  type = "dens_overlay",
  ndraws = 50,
  trim = TRUE,
  re_formula = NULL,
  ...
)
```

## Arguments

- object:

  A `because` fit object.

- resp:

  Character string; the response variable to check. If `NULL`, takes the
  first response.

- type:

  Character string; the type of PPC plot to generate. Supported:
  `"dens_overlay"`, `"hist"`, `"stat"`.

- ndraws:

  Integer; number of posterior draws to use. Defaults to 50.

- re_formula:

  Formula or `NA`; determines which random effects to include in the
  posterior predictions. See `posterior_predict` for details.

- ...:

  Additional arguments passed to `bayesplot` functions.

## Value

A `ggplot` object produced by `bayesplot`.
