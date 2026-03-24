# Initialization Hook

Modules implement this to provide specialized initial values (e.g.
latent states for occupancy).

## Usage

``` r
get_inits_hook(family, data, ...)
```

## Arguments

- family:

  The S3 family list.

- data:

  The model data.

- ...:

  Additional arguments.

## Value

A named list of initial values.
