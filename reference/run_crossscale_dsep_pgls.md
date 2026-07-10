# Run a cross-scale d-sep test via aggregation and PGLS/GLS

Implements the scale-aware d-sep testing strategy for hierarchical
models:

1.  Aggregate the response variable to the predictor's hierarchical
    level (log-mean for Poisson; arithmetic mean for Gaussian).

2.  Select conditioning variables at or above the predictor's level.

3.  Fit a GLS with
    [`ape::corPagel`](https://rdrr.io/pkg/ape/man/corPagel.html)
    (Pagel's lambda estimated by ML) if a phylogenetic tree matching the
    predictor's level is found in `structure`; otherwise fall back to
    ordinary OLS.

4.  Return synthetic MCMC samples (normal approximation to the GLS
    posterior of the focal coefficient) so that the result slots
    directly into because()'s existing downstream summary / plot code.

## Usage

``` r
run_crossscale_dsep_pgls(
  i,
  test_eq,
  cs_info,
  original_data,
  hierarchical_info,
  structure = NULL,
  family = NULL,
  engine = "numpyro",
  n.iter = 1000L,
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

  Named list of structure objects (phylo trees, spatial matrices), as
  passed to
  [`because()`](https://because-pkg.github.io/because/reference/because.md).

- family:

  Named character vector mapping response names to families.

- n.iter:

  Synthetic draws per MCMC chain (default 1000).

- n.chains:

  Number of synthetic chains (default 3).

- quiet:

  Suppress informational messages.

## Value

List with `$samples` (coda::mcmc.list), `$param_map`, `$model`,
`$test_index`.
