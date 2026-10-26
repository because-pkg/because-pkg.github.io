# Get the latent predictor variable name for a family

Extension packages implement this to substitute latent states when a
specialized variable is used as a predictor in downstream equations.

## Usage

``` r
get_latent_predictor_name(family, predictor_name, predictor_index, ...)
```

## Arguments

- family:

  A family object (S3 class `because_family_*`)

- predictor_name:

  Character; the predictor variable name

- predictor_index:

  Character; the indexed predictor string (e.g. `"Species[i]"`)

- ...:

  Additional arguments

## Value

Character; the indexed string to use in the linear predictor
