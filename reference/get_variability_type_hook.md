# Variability Type Hook

Modules implement this to provide default variability metadata (e.g.
'reps' for occupancy).

## Usage

``` r
get_variability_type_hook(family, variable_name, ...)
```

## Arguments

- family:

  The S3 family list.

- variable_name:

  Name of the variable.

- ...:

  Additional arguments.

## Value

Character string ('se', 'reps', or NULL).
