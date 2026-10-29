# Bayesian Homoscedasticity and Residual Dispersion Diagnostics

Assesses variance constancy across fitted predictions (Gelman & Hill
2007) by partitioning fitted values into quantiles and computing
posterior residual standard deviations and variance ratios.

## Usage

``` r
check_homoscedasticity(object, ...)

# S3 method for class 'because'
check_homoscedasticity(object, resp = NULL, n_bins = 5, ndraws = 250, ...)
```

## Arguments

- object:

  A `because` model fit object.

- ...:

  Additional arguments.

- resp:

  Character string; name of the response variable. If `NULL`, takes the
  first response.

- n_bins:

  Integer; number of bins to partition fitted values into. Defaults to 5
  (quintiles).

- ndraws:

  Integer; number of posterior draws to use. Defaults to 250.

## Value

An object of class `"because_homoscedasticity"` containing diagnostic
plots and summary statistics.
