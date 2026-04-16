# Default Method for D-Sep Tree Hook

D-separation independence tests are pure regression tests; they should
not inherit the full phylogenetic or spatial covariance structure from
the main model, which would cause JAGS compilation failures (dimension
mismatches or undefined loop indices). Returning NULL ensures each test
is fitted as a plain mixed-effects regression without structured
covariance.

## Usage

``` r
# Default S3 method
dsep_tree_hook(tree, test_eq, hierarchical_info, levels, ...)
```
