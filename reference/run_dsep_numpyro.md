# Run d-separation tests using NumPyro backend

Evaluates the basis set of conditional independence tests using NumPyro,
handling structure transforms, Python data flattening, cross-scale
tests, and caching.

## Usage

``` r
run_dsep_numpyro(
  equations,
  structures,
  data,
  latent,
  random_terms,
  hierarchical_info,
  all_poly_terms,
  family,
  quiet,
  original_call,
  original_data,
  id_col = NULL,
  variability = NULL,
  adapt_delta = 0.95,
  max_treedepth = 10,
  expand_ordered = FALSE,
  reuse_models = NULL,
  aggregate_crossscale = NULL,
  n.iter = 2000,
  n.burnin = 1000,
  n.chains = 4,
  n.thin = 1,
  parallel = FALSE,
  n.cores = 1,
  dsep_max_obs = 10000,
  prior_scale_fixed = NULL,
  engine = "numpyro",
  structure = NULL,
  because_py = NULL
)
```
