# Bayesian Collinearity and Variance Inflation Diagnostics

Evaluates collinearity among predictors in each equation of a fitted
`because` model using the Posterior Parameter Correlation Matrix and
Leamer's Bayesian Variance Inflation Factor (BVIF).

## Usage

``` r
check_collinearity(object, ...)

# S3 method for class 'because'
check_collinearity(object, equation = NULL, ...)

bvif(object, ...)
```

## Arguments

- object:

  A `because` model fit object.

- ...:

  Additional arguments.

- equation:

  Optional character string or integer specifying a single equation to
  evaluate. If `NULL` (default), evaluates all equations with 2 or more
  predictors.

## Value

An object of class `"because_collinearity"`, containing summary tables
and correlation matrices.
