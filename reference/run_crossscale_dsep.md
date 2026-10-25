# Run a cross-scale d-sep test via coarse-level aggregation

Implements the scale-aware d-sep testing strategy for hierarchical
models:

1.  Aggregate the response variable to the predictor's hierarchical
    level (link-scale mean for non-Gaussian families; arithmetic mean
    for Gaussian).

2.  Aggregate conditioning variables to the predictor's level.

3.  Filter `structure` to any covariance structures matching the
    predictor's level (via `get_order_labels_hook`).

4.  Fit the aggregated single-equation model natively via
    [`because()`](https://because-pkg.github.io/because/reference/because.md)
    and return the posterior samples for the focal coefficient.

## Usage

``` r
run_crossscale_dsep(
  i,
  test_eq,
  cs_info,
  original_data,
  hierarchical_info,
  structure = NULL,
  family = NULL,
  engine = "numpyro",
  n.iter = 1000L,
  n.burnin = 500L,
  n.thin = 1L,
  n.adapt = 500L,
  n.chains = 3L,
  quiet = FALSE
)
```

## Arguments

- i:

  Integer index of the d-sep test.

- test_eq:

  Formula carrying `attr(., "test_var")`.

- cs_info:

  List returned by
  [`detect_crossscale_dsep()`](https://because-pkg.github.io/because/reference/detect_crossscale_dsep.md).

- original_data:

  Hierarchical data list or flat data.frame.

- hierarchical_info:

  List with `$levels`, `$hierarchy`, `$link_vars`.

- structure:

  Named list of structure objects (e.g. spatial matrices, phylogenetic
  trees), as passed to
  [`because()`](https://because-pkg.github.io/because/reference/because.md).

- family:

  Named character vector mapping response names to families.

- engine:

  Bayesian engine (`"jags"`, `"nimble"`, or `"numpyro"`).

- n.iter:

  MCMC iterations per chain (default 1000).

- n.burnin:

  Burn-in iterations per chain (default 500).

- n.thin:

  Thinning interval (default 1).

- n.adapt:

  Adaptation iterations (default 500).

- n.chains:

  Number of MCMC chains (default 3).

- quiet:

  Suppress informational messages.

## Value

List with `$samples` (coda::mcmc.list), `$param_map`, `$model`,
`$test_index`.
