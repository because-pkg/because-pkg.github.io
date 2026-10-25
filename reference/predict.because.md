# Predict Method for because Model Fits

Computes posterior predictions, fitted values, and credible intervals
for a fitted `because` model.

## Usage

``` r
# S3 method for class 'because'
predict(
  object,
  resp = NULL,
  newdata = NULL,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ndraws = NULL,
  re_formula = NULL,
  ...
)
```

## Arguments

- object:

  A `because` fit object.

- resp:

  Character string; name of the response variable to predict. If `NULL`,
  defaults to the first response variable.

- newdata:

  Optional new data frame for prediction. Currently ignored.

- summary:

  Logical; if `TRUE` (default), returns a data frame with summary
  statistics (`Estimate`, `Est.Error`, `Q2.5`, `Q97.5`). If `FALSE`,
  returns the raw draws matrix from
  [`posterior_predict`](https://because-pkg.github.io/because/reference/posterior_predict.md).

- probs:

  Numeric vector of length 2; quantiles for summary intervals. Defaults
  to `c(0.025, 0.975)`.

- ndraws:

  Integer; number of posterior draws to use. Defaults to all draws.

- re_formula:

  Formula or `NA`; determines which random effects to include. See
  [`posterior_predict`](https://because-pkg.github.io/because/reference/posterior_predict.md)
  for details.

- ...:

  Additional arguments passed to
  [`posterior_predict`](https://because-pkg.github.io/because/reference/posterior_predict.md).

## Value

If `summary = TRUE`, a data frame containing:

- Estimate:

  Posterior mean prediction

- Est.Error:

  Posterior standard deviation

- Q2.5:

  Lower credible limit

- Q97.5:

  Upper credible limit

If `summary = FALSE`, an `[ndraws x N_obs]` matrix of posterior draws.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(Y = rnorm(100), X = rnorm(100))
fit <- because(list(Y ~ X), data = df)

# Fitted summary with 95% credible intervals
preds <- predict(fit, resp = "Y")
head(preds)

# Raw posterior draws
draws <- predict(fit, resp = "Y", summary = FALSE)
} # }
```
