# Check if an equation is an auxiliary sub-model equation

Extension packages implement this generic to indicate which equations
should skip standalone likelihood and residual variance generation. Base
`because` still emits linear predictor (`mu_X[i]`) and intercept
(`alpha_X`) priors needed by the sub-model.

## Usage

``` r
is_auxiliary_equation(family, response, all_responses, ...)
```

## Arguments

- family:

  A family object (S3 class `because_family_*`)

- response:

  Character; the response variable name

- all_responses:

  Character vector of all response variable names in the model

- ...:

  Additional arguments

## Value

A list: `list(skip_likelihood = FALSE, skip_variance = FALSE)`
