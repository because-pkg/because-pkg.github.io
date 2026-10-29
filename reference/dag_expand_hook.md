# DAG Expansion Hook

Modules implement this to expand compound latent nodes and define
compound groups for DAG plotting.

## Usage

``` r
dag_expand_hook(family, equations, latent, ...)

# Default S3 method
dag_expand_hook(family, equations, latent, ...)
```

## Arguments

- family:

  The S3 family list or family object.

- equations:

  List of formulas.

- latent:

  Character vector of latent variables.

- ...:

  Additional arguments.

## Value

A list with `equations` (potentially rewired), `latent` (updated
latents), and `compound_groups` (NULL or list of compound groups).
