# Infer Variable Level

Determine which multiscale level a variable belongs to

## Usage

``` r
infer_variable_level(
  var,
  levels,
  data = NULL,
  equations = NULL,
  latent = NULL,
  hierarchy = NULL
)
```

## Arguments

- equations:

  Optional list of model formulas for context

- latent:

  Character vector of latent variables

- hierarchy:

  Optional hierarchy string (e.g. "individual \> obs") to determine
  ordering

## Value

Character, level name
