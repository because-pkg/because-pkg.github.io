# Unified Diagnostic Dashboard for because Models

Generates a comprehensive multi-panel visual diagnostic summary for a
`because` model equation.

## Usage

``` r
check_model(object, ...)

# S3 method for class 'because'
check_model(object, resp = NULL, ndraws = 250, ...)
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

  Integer; number of draws to use for predictive checks. Defaults to
  250.

## Value

A list containing all diagnostic sub-objects and plots, returned
invisibly.
