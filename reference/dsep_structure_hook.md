# D-Sep Structure Hook

Modules implement this to decide if a covariance structure should be
passed to a specific d-separation test.

## Usage

``` r
dsep_structure_hook(structure, test_eq, hierarchical_info, levels, ...)
```

## Arguments

- structure:

  The structure object.

- test_eq:

  The current d-separation test equation.

- hierarchical_info:

  Hierarchical metadata.

- levels:

  Level mapping.

- ...:

  Additional arguments.

## Value

The structure object or NULL.
