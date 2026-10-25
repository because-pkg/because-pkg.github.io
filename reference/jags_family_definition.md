# Define BUGS Family Implementation

Modules implement this method to inject their specific BUGS likelihood
code (for JAGS and NIMBLE) for custom families.

## Usage

``` r
jags_family_definition(family, response, predictors, ...)
```

## Arguments

- family:

  The S3 family object (e.g., because_family_occupancy).

- response:

  The name of the response variable.

- predictors:

  Character vector of predictor names.

- ...:

  Additional arguments.

## Value

A list containing `model_code` (the likelihood block) and
`monitor_params`.
