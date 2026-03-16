# Get default precision prior for a distribution family

Get default precision prior for a distribution family

## Usage

``` r
jags_family_precision_prior(family, param_name, ...)
```

## Arguments

- family:

  A family object

- param_name:

  Name of the parameter (e.g. "tau_e_Abundance")

- ...:

  Additional arguments

## Value

Character string; JAGS prior statement (e.g. "tau_e_Abundance ~
dgamma(1, 1)")
