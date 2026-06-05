# NumPyro Structure Definition

Defines the Python code (JAX logic) to inject for a custom structure
when using `engine = "numpyro"`. Extension packages (e.g.
`because.phybase`) should implement this S3 generic to support their
structures natively in Python.

## Usage

``` r
numpyro_structure_definition(structure, variable_name = "err", ...)
```

## Arguments

- structure:

  The structural object.

- variable_name:

  The name of the variable.

- ...:

  Additional arguments.

## Value

A character string of valid JAX Python code.
