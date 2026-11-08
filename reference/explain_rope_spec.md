# Format human-readable explanations of ROPE intervals

Format human-readable explanations of ROPE intervals

## Usage

``` r
explain_rope_spec(resolved_df, data = NULL)
```

## Arguments

- resolved_df:

  Data frame containing Response, Family, r_min, r_max

- data:

  Optional model data frame to check if Gaussian traits are standardized

## Value

List with 'lines' (character vector) and 'note' (optional string)
