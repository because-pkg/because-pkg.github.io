# Prepare random effects data structures for because models

Expands random effects terms into the data list, creating index vectors,
grouping matrices, and associated loop bounds needed by the model
engines (JAGS, NIMBLE, NumPyro).

## Usage

``` r
prepare_random_effects_data(
  data,
  random_terms,
  equations,
  hierarchical_info,
  is_hierarchical,
  levels,
  family,
  quiet,
  variability = NULL,
  id_col = NULL,
  all_poly_terms = NULL,
  latent = NULL,
  structure = NULL,
  original_data = NULL
)
```
