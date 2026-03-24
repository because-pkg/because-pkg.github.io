# D-Sep Potential Latent Hook

Modules implement this to remove variables from the 'potential latents'
list (e.g. occupancy states).

## Usage

``` r
dsep_potential_latent_hook(family, potential_latents, ...)
```

## Arguments

- family:

  The S3 family list.

- potential_latents:

  Character vector of variables suspected to be latent.

- ...:

  Additional arguments.

## Value

A modified character vector of potential latents.
