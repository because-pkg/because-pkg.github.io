# Needs Zero Inflation Hook

Modules implement this to signal if a variable needs a zero-inflated
likelihood.

## Usage

``` r
needs_zero_inflation_hook(family, variable_name, ...)
```

## Arguments

- family:

  The S3 family list.

- variable_name:

  Name of the response variable.

- ...:

  Additional arguments.

## Value

Logical.
