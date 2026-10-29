# Bayesian Residual Structure and Variance Partition Diagnostics

Evaluates the proportion of variance absorbed by phylogenetic or spatial
random effects in a `because` model via the posterior Variance Partition
Coefficient (VPC).

## Usage

``` r
check_residual_structure(object, ...)

# S3 method for class 'because'
check_residual_structure(object, resp = NULL, ...)
```

## Arguments

- object:

  A `because` model fit object.

- ...:

  Additional arguments.

- resp:

  Character string; name of the response variable. If `NULL`, takes the
  first response.

## Value

An object of class `"because_residual_structure"`.
