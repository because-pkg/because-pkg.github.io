# Default Method for D-Sep Structure Hook

Filters the covariance structure list for a d-separation test sub-model
based on the levels active in that test (or returns `NULL` for flat
models unless overridden by a structure-specific method).

## Usage

``` r
# Default S3 method
dsep_structure_hook(structure, test_eq, hierarchical_info, levels, ...)
```
