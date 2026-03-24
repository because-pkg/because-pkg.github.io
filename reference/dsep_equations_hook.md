# Expanded D-Sep Equation Hook

Modules implement this to add supporting equations (e.g. detection
models) to a d-separation test.

## Usage

``` r
dsep_equations_hook(family, equations, dsep_equations, ...)
```

## Arguments

- family:

  The S3 family list.

- equations:

  The full list of model equations.

- dsep_equations:

  The current list of equations for the test.

- ...:

  Additional arguments.

## Value

A modified list of equations.
