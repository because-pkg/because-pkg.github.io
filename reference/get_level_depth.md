# Get Depth of a Level in Hierarchy

Correctly handles multi-chain hierarchies (semicolons). Returns the
1-indexed depth (1 = coarsest). If a level appears in multiple chains,
returns the maximum depth.

## Usage

``` r
get_level_depth(lvl, hierarchy_str)
```

## Arguments

- lvl:

  Level name to find

- hierarchy_str:

  Hierarchy string

## Value

Integer depth
