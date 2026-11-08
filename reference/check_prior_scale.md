# Check for unstandardized continuous variables with large variance

Inspects continuous (Gaussian) response variables in the model
equations. If any have large empirical standard deviations (SD \> 20)
and no custom prior was provided for that variable, issues an
informative warning alerting the user that default priors assume
standardized variables.

## Usage

``` r
check_prior_scale(data, equations, family = NULL, priors = NULL, quiet = FALSE)
```

## Arguments

- data:

  Model data (data.frame or list)

- equations:

  List of model equations (formulas)

- family:

  Named vector or list of distribution families

- priors:

  List of custom priors provided by the user (or NULL)

- quiet:

  Logical; if TRUE, suppress warnings
