# Define BUGS Structure Implementation

Modules implement this method to inject their specific BUGS code (for
JAGS and NIMBLE) for covariance structures.

## Usage

``` r
jags_structure_definition(structure, variable_name = "err", ...)
```

## Arguments

- structure:

  The structure object (e.g., phylo, list, custom class).

- variable_name:

  Name of the error variable (default "err").

- ...:

  Additional arguments.

## Value

A list containing `setup_code` (e.g. priors, matrix inversion) and
`error_prior` (the likelihood/prior definition for residuals).
