# Parse Hierarchy from Random Effects

Extract multiscale nesting structure from random effects formula

## Usage

``` r
parse_hierarchy_from_random(random, data = NULL)
```

## Arguments

- random:

  Formula specifying random effects

- data:

  List of data.frames at different hierarchical levels (optional,
  currently unused)

## Value

Character string hierarchy (e.g., "site_year \> individual") or NULL
