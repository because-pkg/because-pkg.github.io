# Bayesian Normality Check and Q-Q Envelope

Evaluates distributional and normality assumptions using Randomized
Quantile Residuals (Dunn & Smyth 1996) and Posterior Predictive Credible
Envelopes (Landwehr et al. 1984).

## Usage

``` r
check_normality(object, ...)

# S3 method for class 'because'
check_normality(object, resp = NULL, ndraws = 250, ...)

plot_qq(object, ...)

# S3 method for class 'because'
plot_qq(object, ...)
```

## Arguments

- object:

  A `because` model fit object.

- ...:

  Additional arguments.

- resp:

  Character string; name of the response variable. If `NULL`, takes the
  first response.

- ndraws:

  Integer; number of posterior predictive draws to simulate for the
  envelope. Defaults to 250.

## Value

An object of class `"because_normality"` (for `check_normality`) or a
`ggplot` object (for `plot_qq`).
