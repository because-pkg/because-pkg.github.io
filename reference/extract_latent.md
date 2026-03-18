# Extract Latent Parameters from a because Model

This function extracts posterior distributions of latent variables such
as occupancy probability (psi), detection probability (p), and latent
occupancy state (z). It is designed to handle multispecies models by
mapping stacked results back to species and sites.

## Usage

``` r
extract_latent(object, type = "occupancy", variables = NULL)
```

## Arguments

- object:

  A `because` model object.

- type:

  Character; the type of latent variable to extract. Currently supports
  "occupancy".

- variables:

  Optional character vector specifying which latent variables to extract
  (e.g., `c("psi", "p", "z")`). If `NULL`, extracts all available for
  the given type.

## Value

A tidy data frame (tibble) containing:

- `Variable`: The base name of the latent parameter (e.g., "psi").

- `Response`: The response variable name in the model (e.g., "Y").

- `SpeciesID`: The species identifier (for multispecies models).

- `SiteID`: The site identifier.

- `Mean, SD, Q2.5, Q50, Q97.5`: Posterior summary statistics.

## Details

For occupancy models, the function looks for parameters starting with
`psi_`, `p_`, and `z_`. It uses the internal `results$data` or
`results$species_order` to map indices back to meaningful identifiers.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming a multispecies occupancy model was fit:
fit <- because(list(Y ~ 1), data = data, family = list(Y = "occupancy"))
latents <- extract_latent(fit, type = "occupancy", variables = c("psi", "p"))
head(latents)
} # }
```
