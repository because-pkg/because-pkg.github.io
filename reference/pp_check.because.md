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

- trim:

  Logical; if TRUE, zooms plot to observed data range.

- re_formula:

  Formula or `NA`; determines which random effects to include in the
  posterior predictions. See
  [`posterior_predict`](https://because-pkg.github.io/because/reference/posterior_predict.md)
  for details.

- ...:

  Additional arguments passed to `bayesplot` functions.

## Value

A `ggplot` object produced by `bayesplot`.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(Y = rnorm(100), X = rnorm(100))
fit <- because(list(Y ~ X), data = df)

# Default density overlay check
pp_check(fit)

# Histogram check for a specific response
pp_check(fit, resp = "Y", type = "hist", ndraws = 30)

# Test statistic check (e.g., mean)
pp_check(fit, type = "stat", stat = "mean")
} # }
```
