# Aggregate Multiscale Data to a Target Resolution

Collapses finer-level variables to a coarser grain by taking the mean.
Used for cross-scale d-separation tests where tests must be locked to
the resolution of the coarsest variable involved.

## Usage

``` r
aggregate_multiscale_data(data, levels, hierarchy, target_lvl)
```

## Arguments

- data:

  List of dataframes

- levels:

  List mapping variables to levels

- hierarchy:

  Hierarchy string (e.g. "Year \> Individual")

- target_lvl:

  The resolution to aggregate everything to

## Value

List of dataframes truncated/aggregated to stop at target_lvl
