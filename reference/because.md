# Run a Bayesian Structural Equation Model

This function fits a Bayesian ...

## Usage

``` r
because(
  equations,
  data,
  id_col = NULL,
  structure = NULL,
  tree = NULL,
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
  reuse_models = NULL
)
```

## Arguments

- equations:

  A list of model formulas describing the structural equation model.

- data:

  A data.frame containining the variables in the model. If using
  hierarchical models (see hierarchical section below), this can also be
  a list of data frames.

- id_col:

  Character string specifying the column name in a data.frame containing
  unit identifiers (e.g., individuals, sites, or species names). This is
  used to:

  - Match data rows to external structure labels (e.g., tip labels in
    phylogenetic trees).

  - Link data to external spatial or custom covariance matrices.

  **Note**: For standard random effects models (e.g.,
  `random = ~(1|group)`) where no external structure is provided, this
  argument is **not required**. The grouping column is read directly
  from the data.

  If `NULL` (default): uses meaningful row names if available. Ignored
  when `data` is already a list.

- structure:

  The covariance structure for the model. Accepts:

  - `NULL`: Independent model (Standard SEM, no covariance structure).

  - `matrix`: Custom covariance or precision matrix (e.g., spatial
    connectivity, kinship).

  - Custom objects: Supported via extension packages (e.g., phylogenetic
    trees from because.phybase or spatial structures).

- tree:

  (Deprecated alias for `structure`). Use `structure` instead for new
  code.

- engine:

  Character string specifying the inference engine to use. Supported
  values:

  - `"jags"` (default): Use Just Another Gibbs Sampler (via rjags).

  - `"nimble"`: Use the compiled C++ backend (via NIMBLE). Offers
    significant speedups for complex models, parallel execution, and
    marginalized likelihoods.

- monitor:

  Parameter monitoring mode. Options:

  - `"interpretable"` (default): Monitor only scientifically meaningful
    parameters: intercepts (alpha), regression coefficients (beta),
    phylogenetic signals (lambda) for responses, and WAIC terms.
    Excludes variance components (tau) and auxiliary predictor
    parameters.

  - `"all"`: Monitor all model parameters including variance components
    and implicit equation parameters.

  - Custom vector: Provide a character vector of specific parameter
    names to monitor.

  - `NULL`: Auto-detect based on model structure (equivalent to
    "interpretable").

- nimble_samplers:

  (NIMBLE-only) A named list specifying custom samplers for specific
  model nodes. Example: `nimble_samplers = list(beta_X_Y = "slice")`.
  Common sampler types include:

  - `"RW"`: Scalar Random-Walk Metropolis-Hastings.

  - `"RW_block"`: Multivariate Random-Walk Metropolis-Hastings.

  - `"slice"`: Scalar slice sampler.

  - `"AF_slice"`: Automated Factor Slice Sampler (multivariate slice).

  - `"categorical"`: Specialized discrete sampler for
    Multinomial/Ordinal choices.

- n.chains:

  Number of MCMC chains (default = 3).

- n.iter:

  Total number of MCMC iterations (default = 12500).

- n.burnin:

  Number of burn-in iterations (default = n.iter / 5).

- n.thin:

  Thinning rate (default = 10).

- DIC:

  Logical; whether to compute DIC using `dic.samples()` (default =
  TRUE). **Note**: DIC penalty will be inflated for models with
  measurement error or repeated measures because latent variables are
  counted as parameters (penalty ~ structural parameters + N). For model
  comparison, use WAIC or compare mean deviance across models with
  similar structure.

- WAIC:

  Logical; whether to sample values for WAIC and deviance (default =
  FALSE). WAIC is generally more appropriate than DIC for hierarchical
  models with latent variables.

- n.adapt:

  Number of adaptation iterations (default = n.iter / 5).

- quiet:

  Logical; suppress JAGS output (default = FALSE).

- verbose:

  Logical; if `TRUE`, print generated JAGS model code and data names
  (default = FALSE).

- dsep:

  Logical; if `TRUE`, evaluate the model's global fit using d-separation
  (basis set) path analysis. This identifies the complete set of
  independence claims implied by the DAG and tests each one via Bayesian
  inference. Results can be explored via
  [`summary()`](https://rdrr.io/r/base/summary.html) or visualized using
  [`plot_dsep()`](https://because-pkg.github.io/because/reference/plot_dsep.md).

- variability:

  Optional specification for variables with measurement error or
  within-species variability. **Global Setting**:

  - `"reps"`: Applies repeat-measures modeling to **all** continuous
    variables in the equations (except grouping variables). Expects
    `X_obs` matrix or long-format data.

  - `"se"`: Applies measurement error modeling to **all** continuous
    variables. Expects `X_se` columns.

  **Manual Specification** (Named Vector/List):

  - Simple: `c(X = "se", Y = "reps")` - mixed types

  - Custom columns: `list(X = list(type = "se", se_col = "X_sd"))`

  - For SE: `se_col` (SE column), `mean_col` (mean column, optional)

  - For reps: `obs_col` (observations matrix column)

  **Auto-Detection**: If not specified, the package attempts to detect
  variability based on column names:

  - `X_se` -\> type="se"

  - `X_obs` or matrix column -\> type="reps"

- family:

  Optional named character vector specifying the family/distribution for
  response variables. Additional families (e.g., `"occupancy"`) are
  provided by the because.detection extension package. Example:
  `family = c(Gregarious = "binomial")`.

- distribution:

  Deprecated alias for `family`.

- latent:

  Optional character vector of latent (unmeasured) variable names. If
  specified, the model will account for induced correlations among
  observed variables that share these latent common causes.

- latent_method:

  Method for handling latent variables (default = "correlations").

  - `"correlations"`: MAG approach - marginalize latent variables and
    estimate induced correlations (`rho`) between observed variables
    that share latent parents.

  - `"explicit"`: Model latent variables as JAGS nodes and estimate
    structural paths from latents to observed variables.

- standardize_latent:

  Logical; if `TRUE` and `latent_method = "explicit"`, adds standardized
  priors (`N(0,1)`) to latent variables to identify scale and location.
  This improves convergence and makes regression coefficients
  interpretable as standardized effects. Only applicable when using
  explicit latent variable modeling (default = TRUE).

- parallel:

  Logical; if `TRUE`, run MCMC chains in parallel (default = FALSE). For
  standard SEM (`dsep = FALSE`), this runs the chains on different
  cores. For d-separation testing (`dsep = TRUE`), this instead runs
  individual independence tests on different cores to maximize
  throughput. Each sub-test runs its chains sequentially. Note: Requires
  `n.cores > 1` to take effect.

- n.cores:

  Integer; number of CPU cores to use for parallel chains (default = 1).
  Only used when `parallel = TRUE`.

- cl:

  Optional cluster object for parallel execution.

- ic_recompile:

  Logical; if `TRUE` and `parallel = TRUE`, recompile the model after
  parallel chains to compute DIC/WAIC (default = TRUE). This adds a
  small sequential overhead but enables information criteria
  calculation.

- random:

  Optional formula or list of formulas specifying global random effects
  applied to all equations (e.g. `~(1|species)`).

- levels:

  (Hierarchical Data) A named list mapping variables to their hierarchy
  levels. Required if `data` is a list of data frames (hierarchical
  format). Example: `list(individual = c("y", "x"), site = c("z"))`.

- hierarchy:

  (Hierarchical Data) Character string describing the topological
  ordering of levels (e.g., `"site > individual"`). Required for
  hierarchical data if not fully inferred from random effects.

- link_vars:

  (Hierarchical Data) Optional named character vector specifying
  variables used to link data levels (e.g. `c(site = "site_id")`).

- fix_residual_variance:

  Optional named vector for fixing residual variances. Useful for
  handling non-identified models or specific theoretical constraints.
  Example: `c(response_var = 1)`.

- priors:

  Optional named list of character strings specifying custom priors for
  specific parameters. By default, `because` uses "boundary-avoiding"
  regularizing priors for hierarchical variance components:

  - **Gaussian**: `dgamma(1, 1)` (vague) for residual and random effect
    precisions.

  - **Non-Gaussian**: `dgamma(10, 10)` for overdispersion and
    hierarchical random effects (\\\tau\\). This regularizing prior
    (Gelman 2006, McElreath 2020) prevents numerical overflow in models
    with exponential or logit link functions by constraining the sampler
    away from astronomically large variances during adaptation.

  Example:
  `list(alpha_Response = "dnorm(0, 0.001)", tau_e_Response = "dgamma(1, 1)")`.

- reuse_models:

  List of previously fitted 'because' models to scan for reusable
  d-separation test results. If a test in the current run matches a test
  in a reused model (same formula), the result is copied instead of
  re-running JAGS. **Note**: Ensuring that the data is consistent is the
  user's responsibility.

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
fit_med <- because(equations, data = match_df)
} # }
```
