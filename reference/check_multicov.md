# Check Multi-Structure Covariance Convergence and Sensitivity Diagnostics

Evaluates structure sampling, within-structure MCMC convergence,
within-group structure sensitivity, and between-group divergence
(`Rhat_group` / `Rhat_backbone`), alongside an exact Law of Total
Variance decomposition partitioning posterior parameter uncertainty
(total \\\text{SD}^2\\) into within-tree estimation uncertainty
(unexplained biological variance / sample size), within-group structure
variation, and between-group hypothesis variation for any `because`
model fitted with a collection or ensemble of covariance structures
(e.g., phylogenetic trees/backbones, spatial distance kernels, genomic
relationship matrices, or network graphs).

## Usage

``` r
check_multicov(model, plot = TRUE, show_internal = FALSE, ...)
```

## Arguments

- model:

  A `because` model object fitted with multiple candidate covariance
  structures.

- plot:

  Logical. If `TRUE` (default), displays a publication-ready 2-panel
  diagnostic figure featuring hypothesis sensitivity forest plot (left)
  and 100% stacked uncertainty variance decomposition (right).

- show_internal:

  Logical. If `FALSE` (default), suppresses internal additive simulation
  noise parameters (e.g. `sigma_Brain`, `sigma_shared_*`) for pure
  latent indicator variables in marginalized latent (MAG) models.

- ...:

  Additional arguments passed to plotting functions or metadata
  overrides (e.g., `group = c(...)` to override group/backbone labels).

## Value

An invisible object of class `check_multicov` containing:

- param_diagnostics:

  Data frame of parameter sensitivity \\\widehat{R}\\ across candidate
  structures/trees, groups/backbones, and overall ensemble omnibus
  \\\widehat{R}\_{\text{all}}\\.

- variance_decomposition:

  Exact Law of Total Variance table decomposing total parameter
  posterior variance (\\\text{SD}^2\\) into within-structure estimation
  uncertainty, between-structure variation, and between-group hypothesis
  variation.

- rhat_per_tree:

  Matrix of within-structure MCMC \\\widehat{R}\\ values for each
  parameter on each candidate tree/structure. Kept unprinted by default
  to avoid console clutter, but directly indexable on demand.

- structure_weights_pct:

  Posterior weights per candidate structure/tree.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit with multiple candidate covariance structures (e.g. multiPhylo ensemble)
fit <- because(list(y ~ x), data = my_data, structure = my_cov_list, engine = "numpyro")
chk <- check_multicov(fit)

# Inspect within-tree MCMC convergence for specific parameters or trees on demand:
chk$rhat_per_tree["beta_y_x", ]   # Within-tree Rhat across all candidate trees for beta_y_x
chk$rhat_per_tree[, "Tree_1"]     # Within-tree Rhat for all parameters on Tree 1
} # }
```
