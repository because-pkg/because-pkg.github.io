# NIMBLE Family Optimization

Modules implement this to provide specialized NIMBLE optimizations (e.g.
marginalization).

## Usage

``` r
nimble_family_optimization(family, model_string, ...)
```

## Arguments

- family:

  The S3 family object.

- model_string:

  The current model code string.

- ...:

  Additional arguments.

## Value

A list with `model_string` and `nimble_functions` (list of
nimbleFunction).
