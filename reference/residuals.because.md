# Extract Fitted Values and Residuals for because Models

Bayesian methods for extracting expected fitted values and residuals
from fitted `because` model objects.

## Usage

``` r
# S3 method for class 'because'
fitted(
  object,
  resp = NULL,
  newdata = NULL,
  ndraws = NULL,
  re_formula = NULL,
  summary = TRUE,
  ...
)

# S3 method for class 'because'
residuals(
  object,
  resp = NULL,
  type = c("quantile", "response", "pearson"),
  summary = TRUE,
  ndraws = 250,
  re_formula = NULL,
  ...
)
```

## Arguments

- object:

  A `because` model fit object.

- resp:

  Character string; name of the response variable. If `NULL`, takes the
  first response.

- newdata:

  Optional new dataset (currently evaluated within-sample).

- ndraws:

  Integer; number of posterior draws to use. Defaults to all draws (or
  250 for quantile residuals).

- re_formula:

  Formula or `NA`; determines whether random/phylogenetic effects are
  included. If `NULL` (default), conditional on random effects. If `NA`,
  marginal prediction without random effects.

- summary:

  Logical; if `TRUE` (default), returns a summarized vector (posterior
  mean or median). If `FALSE`, returns the full `[ndraws x N_obs]`
  posterior matrix.

- type:

  Character; type of residual: `"quantile"` (default), `"response"`, or
  `"pearson"`.

- ...:

  Additional arguments.

## Value

A numeric vector (if `summary = TRUE`) or matrix of posterior draws (if
`summary = FALSE`).

## References

Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, 5(3), 236-244.

Hartig, F. (2022). DHARMa: Residual Diagnostics for Hierarchical
(Multi-Level / Mixed) Regression Models (R package version 0.4.6).
r-project.org

Gelman, A., & Hill, J. (2007). *Data Analysis Using Regression and
Multilevel/Hierarchical Models*. Cambridge University Press.
