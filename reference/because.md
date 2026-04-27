# Run a Bayesian Structural Equation Model (Because)

Fits a Bayesian Structural Equation Model (SEM) using JAGS or NIMBLE.
Supports multi-level (hierarchical) data, custom covariance structures
(phylogenetic, spatial, etc.), missing data imputation, and d-separation
global fit testing.

## Usage

``` r
because(
  equations,
  data,
  id_col = NULL,
  structure = NULL,
  engine = "jags",
  monitor = "interpretable",
  nimble_samplers = NULL,
  n.chains = 3,
  n.iter = 12500,
  n.burnin = floor(n.iter/5),
  n.thin = 10,
  DIC = TRUE,
  WAIC = FALSE,
  n.adapt = floor(n.iter/5),
  quiet = FALSE,
  verbose = FALSE,
  dsep = FALSE,
  variability = NULL,
  family = NULL,
  distribution = NULL,
  latent = NULL,
  latent_method = "correlations",
  standardize_latent = TRUE,
  fix_latent = "loading",
  parallel = FALSE,
  n.cores = parallel::detectCores() - 1,
  cl = NULL,
  ic_recompile = FALSE,
  random = NULL,
  levels = NULL,
  hierarchy = NULL,
  link_vars = NULL,
  fix_residual_variance = NULL,
  priors = NULL,
  reuse_models = NULL,
  expand_ordered = FALSE,
  structure_multi = NULL,
  structure_levels = NULL,
  ...
)
```

## Arguments

- equations:

  A list of model formulas describing the structural equation model.

- data:

  A data.frame containing the variables in the model. If using
  hierarchical models, this can also be a list of data frames.

- id_col:

  Character string specifying the column name containing unit
  identifiers (e.g., species names). If NULL, uses row names.

- structure:

  Optional structural object (e.g., matrix, tree) for correlated
  residuals.

- engine:

  Bayesian engine to use: "jags" (default) or "nimble".

- monitor:

  Character; "interpretable" (default) or "all" parameters to monitor.

- nimble_samplers:

  Optional named list of user-specified NIMBLE samplers.

- n.chains:

  Integer; number of MCMC chains (default = 3).

- n.iter:

  Integer; total iterations per chain (default = 12500).

- n.burnin:

  Integer; burn-in iterations (default = 20% of n.iter).

- n.thin:

  Integer; thinning interval (default = 10).

- DIC:

  Logical; compute DIC? (default = TRUE).

- WAIC:

  Logical; compute WAIC? (default = FALSE).

- n.adapt:

  Integer; adaptation iterations (default = 20% of n.iter).

- quiet:

  Logical; suppress MCMC progress messages? (default = FALSE).

- verbose:

  Logical; print verbose debug information? (default = FALSE).

- dsep:

  Logical; perform d-separation independence tests? (default = FALSE).

- variability:

  Optional specification for variables with measurement error or
  within-species variability.

- family:

  Optional named character vector specifying the family/distribution for
  response variables.

- distribution:

  Deprecated alias for `family`.

- latent:

  Optional character vector of latent (unmeasured) variable names.

- latent_method:

  Method for handling latent variables ("correlations" or "explicit").

- standardize_latent:

  Logical; standardize latents to unit variance?

- fix_latent:

  Identification strategy for latent variables ("loading" or "sign").

- parallel:

  Logical; run chains in parallel? (default = FALSE).

- n.cores:

  Integer; number of CPU cores to use.

- cl:

  Optional cluster object for parallel execution.

- ic_recompile:

  Logical; recompile model after parallel chains for IC?

- random:

  Optional formula for global random effects (e.g. ~(1\|species)).

- levels:

  (Hierarchical) Named list mapping variables to hierarchy levels.

- hierarchy:

  (Hierarchical) Topological ordering of levels (e.g., "site \>
  individual").

- link_vars:

  (Hierarchical) Named vector specifying variables used to link levels.

- fix_residual_variance:

  Optional named vector for fixing residual variances.

- priors:

  Optional named list of custom priors.

- reuse_models:

  List of previously fitted 'because' models to scan for reusable
  results.

- expand_ordered:

  Logical; expand ordered factors into polynomial contrasts?

- structure_multi:

  List of multiple structures for uncertainty.

- structure_levels:

  Mapping of structures to hierarchy levels.

- ...:

  Additional arguments passed to because_model.

## Value

An object of class `"because"` containing:

- samples:

  MCMC samples (mcmc.list).

- parameter_map:

  Data frame mapping parameter names to model variables.

- model:

  JAGS model code.

- summary:

  Summary of posterior samples.

- dsep_results:

  List of d-separation test results (if dsep=TRUE).

- DIC:

  Deviance Information Criterion (if DIC=TRUE).

- WAIC:

  Watanabe-Akaike Information Criterion (if WAIC=TRUE).

## Examples

``` r
if (FALSE) { # \dontrun{
# A simple linear regression
df <- data.frame(Y = rnorm(100), X = rnorm(100))
fit <- because(list(Y ~ X), data = df)

# A mediation model
equations <- list(
  M ~ X,
  Y ~ M + X
)
fit_med <- because(equations, data = df)
} # }
```
