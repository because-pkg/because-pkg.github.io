# Get Data for Variables

Determine the finest grain level needed for a set of variables and
return the appropriate dataset (with joining if needed)

## Usage

``` r
get_data_for_variables(
  variables,
  data,
  levels,
  hierarchy,
  link_vars,
  equations = NULL,
  latent = NULL
)
```

## Arguments

- variables:

  Character vector of variable names

- data:

  List of data.frames at different levels

- levels:

  List mapping variable names to level names

- hierarchy:

  Character string specifying nesting

- link_vars:

  Character vector of linking variables

## Value

data.frame
