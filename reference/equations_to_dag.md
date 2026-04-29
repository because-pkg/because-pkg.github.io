# Convert because equations to dagitty-style adjacency matrix

When `deterministic_terms` is supplied (a list returned by
`extract_deterministic_terms`), interaction and
[`I()`](https://rdrr.io/r/base/AsIs.html) terms are kept as **explicit
intermediate nodes** in the DAG rather than collapsed to their component
variables. This is required to produce the correct conditional
independence basis set following Geiger, Verma & Pearl (1990), which
extends d-separation to handle deterministic nodes.

## Usage

``` r
equations_to_dag(equations, exclude_vars = NULL, deterministic_terms = NULL)
```

## Arguments

- equations:

  List of formulas

- exclude_vars:

  Character vector of variable names to exclude (e.g., grouping
  variables)

- deterministic_terms:

  Optional named list returned by `extract_deterministic_terms`. Each
  element must have `$original` (the R term string, e.g. `"BM:M"`) and
  `$internal_name` (the JAGS-safe node name, e.g. `"BM_x_M"`).

## Value

Named adjacency matrix in binary format

## References

Geiger, D., Verma, T., & Pearl, J. (1990). Identifying independence in
Bayesian Networks. *Networks*, 20(5), 507–534.
