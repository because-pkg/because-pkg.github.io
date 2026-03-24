# Normalize Equations Hook

Modules implement this to normalize specialized aliases in formulas
(e.g. psi_Species -\> Species).

## Usage

``` r
normalize_equations_hook(family, equations, ...)
```

## Arguments

- family:

  The S3 family list.

- equations:

  The list of structural formulas.

- ...:

  Additional arguments.

## Value

A modified list of formulas.
