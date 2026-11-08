# Check Multi-Structure Covariance Convergence and Sensitivity Diagnostics

Evaluates structure sampling, within-structure MCMC convergence
(`Rhat_MCMC`), within-group structure sensitivity, and between-group
divergence (`Rhat_group`), alongside an exact Law of Total Variance
decomposition partitioning posterior uncertainty into MCMC estimation
variance, within-group structure variation, and between-group conflict
for any `because` model fitted with a collection or ensemble of
covariance structures (e.g., phylogenetic trees/backbones, spatial
distance kernels, genomic relationship matrices, or network graphs).

## Usage

``` r
check_multicov(model, plot = TRUE, ...)
```

## Arguments

- model:

  A `because` model object fitted with multiple candidate covariance
  structures.

- plot:

  Logical. If `TRUE` (default), displays a publication-ready 2-panel
  diagnostic figure featuring convergence/sensitivity forest plot (left)
  and 100% stacked uncertainty variance decomposition (right).

- ...:

  Additional arguments passed to plotting functions or metadata
  overrides (e.g., `group = c(...)` to override group/backbone labels).

## Value

An invisible object of class `check_multicov` containing per-chain
structure statistics, structure frequencies, parameter-level split
`Rhat_MCMC` vs. `Rhat_group` diagnostics, and the exact uncertainty
`variance_decomposition` table.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit with multiple candidate covariance structures
fit <- because(list(y ~ x), data = my_data, structure = my_cov_list, engine = "numpyro")
check_multicov(fit)
} # }
```
