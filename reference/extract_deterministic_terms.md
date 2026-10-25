# Extract Deterministic Terms from Formulas

Scans a list of formulas for terms that require deterministic node
definitions (interactions, [`I()`](https://rdrr.io/r/base/AsIs.html)
expressions, etc.) across model building, d-separation, and DAG
construction.

## Usage

``` r
extract_deterministic_terms(equations)
```

## Arguments

- equations:

  List of formulas

## Value

List of deterministic term definitions
