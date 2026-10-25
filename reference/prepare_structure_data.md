# Prepare Structure Data

Modules implement this to process the structure object into data for the
model engine (JAGS, NIMBLE, or NumPyro).

## Usage

``` r
prepare_structure_data(structure, data, ...)
```

## Arguments

- structure:

  The structure object.

- data:

  The model data.

- ...:

  Additional arguments.

## Value

A named list of data to be passed to the model engine (e.g.,
`list(data_list = list(Prec = ...))`).
