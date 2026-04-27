# Check if a d-sep test is a cross-hierarchy test (trivially satisfied)

A test Response *\|\|* FocalPredictor \| \\..\\ is "cross-hierarchy"
when the response and the focal predictor live in orthogonal
hierarchical branches (e.g., one is species-level, the other is
site-level) AND the conditioning set does not contain a variable from a
level that connects the two branches (typically the observation level).
Such tests are trivially independent by design of the hierarchical data
structure and cannot be run in JAGS because no cross-level index exists
between orthogonal branches.

## Usage

``` r
is_cross_hierarchy_test(test_eq, hierarchical_info)
```

## Arguments

- test_eq:

  A formula with attribute "test_var" naming the focal predictor.

- hierarchical_info:

  List with \$levels (named list) and \$hierarchy (string).

## Value

TRUE if the test should be skipped; FALSE otherwise.
