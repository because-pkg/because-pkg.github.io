# Default Method for D-Sep Tree Hook

Filters the covariance structure list for a d-separation test sub-model
based on the levels active in that test (or returns `NULL` for flat
models unless overridden by a structure-specific method).

## Usage

``` r
# Default S3 method
dsep_tree_hook(tree, test_eq, hierarchical_info, levels, ...)
```
