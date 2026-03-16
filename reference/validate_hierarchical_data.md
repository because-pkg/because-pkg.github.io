# Validate Hierarchical Data Structure

Validate Hierarchical Data Structure

## Usage

``` r
validate_hierarchical_data(
  data,
  levels,
  hierarchy,
  link_vars,
  latent_vars = NULL
)
```

## Arguments

- data:

  List of data.frames at different hierarchical levels

- levels:

  List mapping variable names to level names

- hierarchy:

  Character string specifying nesting (e.g., "site_year \> individual")

- link_vars:

  Character vector of variables that link levels
