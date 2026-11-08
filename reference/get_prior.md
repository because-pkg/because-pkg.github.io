# Get prior for a parameter (custom override or default)

Get prior for a parameter (custom override or default)

## Usage

``` r
get_prior(param_name, type = "beta", default = NULL, priors = NULL)
```

## Arguments

- param_name:

  Name of the parameter

- type:

  Type of parameter ("beta", "alpha", etc.)

- default:

  Default prior string

- priors:

  Named list of custom priors

## Value

Formatted BUGS prior statement
