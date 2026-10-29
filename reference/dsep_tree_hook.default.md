# Default Method for D-Sep Tree Hook

Filters the covariance structure list for a d-separation test sub-model
based on the levels active in that test (or returns `NULL` for flat
models unless overridden by a structure-specific method).

## Usage

``` r
# Default S3 method
dsep_tree_hook(structure, test_eq, hierarchical_info, levels, ...)
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
