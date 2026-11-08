# Run d-separation tests using JAGS or NIMBLE backends

Evaluates the basis set of conditional independence tests using JAGS or
NIMBLE, supporting sequential or parallel execution, cross-scale tests,
model reuse, and combining resulting MCMC samples.

## Usage

``` r
run_dsep_jags_nimble(
  engine,
  equations,
  data,
  family,
  structures,
  structure_obj,
  hierarchical_info,
  is_hierarchical,
  random_terms,
  levels,
  multiscale,
  link_vars,
  latent,
  latent_method,
  id_col,
  variability,
  all_poly_terms,
  fixed_equations_temp,
  induced_cors_in = NULL,
  dsep_max_obs,
  aggregate_crossscale,
  parallel,
  n.cores,
  n.chains,
  n.iter,
  n.burnin,
  n.thin,
  n.adapt,
  ic_recompile,
  fix_residual_variance,
  quiet,
  priors,
  monitor_mode,
  expand_ordered,
  nimble_samplers,
  adapt_delta,
  max_treedepth,
  prior_scale_fixed,
  verbose,
  family_obj = NULL,
  hierarchy = NULL,
  original_call = NULL,
  original_data = NULL,
  random = NULL,
  response_vars_with_na = NULL,
  reuse_models = FALSE,
  structure = NULL,
  WAIC = FALSE,
  DIC = FALSE
)
```
