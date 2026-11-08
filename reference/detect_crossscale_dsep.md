# Detect whether a d-sep test crosses hierarchical scales in the same chain

A test is "cross-scale" when the focal predictor is at a coarser level
than the response within the SAME hierarchical chain (not orthogonal
branches). E.g., a site- or species-level variable predicting an
obs-level response in "site \> survey \> obs; species \> obs".

## Usage

``` r
detect_crossscale_dsep(test_eq, hierarchical_info)
```

## Arguments

- test_eq:

  A formula carrying the attribute `test_var`.

- hierarchical_info:

  List with `$levels`, `$hierarchy`, `$link_vars`.

## Value

A list. If cross-scale: `$is_crossscale = TRUE` plus `$response`,
`$response_level`, `$predictor_level`, `$test_var`. Otherwise
`$is_crossscale = FALSE`.
