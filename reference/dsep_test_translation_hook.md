# D-Sep Test Translation Hook

Modules implement this to translate test equations (e.g. psi_Species ~ X
to Species ~ X).

## Usage

``` r
dsep_test_translation_hook(family, test_eq, ...)
```

## Arguments

- family:

  The S3 family list.

- test_eq:

  The d-separation test equation.

- ...:

  Additional arguments.

## Value

A modified test equation.
