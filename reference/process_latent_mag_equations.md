# Process Latent Variable Structural Equations for MAG Logic

Marginalizes latent variables from structural equations, computes
induced correlations if needed, and injects intercept-only models for
disconnected variables.

## Usage

``` r
process_latent_mag_equations(
  latent,
  latent_method,
  equations,
  random_terms,
  hierarchical_info,
  family,
  quiet = FALSE,
  induced_cors = NULL,
  dsep = FALSE
)
```
