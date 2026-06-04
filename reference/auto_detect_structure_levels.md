# Auto-Detect Structure Levels

Maps covariance structures (e.g. phylo, spatial) to hierarchical levels
by matching dimensions and ID values.

## Usage

``` r
auto_detect_structure_levels(structure, hierarchical_info, quiet = FALSE)
```

## Arguments

- structure:

  List of structure objects

- hierarchical_info:

  List containing 'data', 'levels', etc.

- quiet:

  Logical

## Value

Named list mapping structure names to level names
