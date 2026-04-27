# Define JAGS Family Implementation

Define JAGS Family Implementation

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
