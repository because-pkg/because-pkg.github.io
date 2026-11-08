# Default and Resolved Region of Practical Equivalence (ROPE)

Internal helpers for resolving family- and variable-specific ROPE
intervals, formatting human-readable explanations, and providing
principled defaults.

## Usage

``` r
get_family_default_rope(fam)
```

## Arguments

- fam:

  Character; family name (e.g. "gaussian", "poisson", "binomial")

## Value

Numeric vector of length 2: c(min, max)
