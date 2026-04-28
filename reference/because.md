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

  A list of R formulas describing the structural model. Each formula
  represents a causal path (e.g., `Y ~ X1 + X2`). Categorical predictors
  are automatically handled.

- data:

  A data frame containing all variables. For **hierarchical models**,
  this must be a named list of data frames (one per level).

- id_col:

  Character string for the column identifying units (e.g. species
  names).

- structure:

  Optional covariance structure (e.g., a phylogenetic tree or spatial
  matrix) for modeling correlated residuals.

- engine:

  Bayesian backend: `"jags"` (default) or `"nimble"`.

- n.chains:

  Number of independent MCMC chains (default = 3).

- n.iter:

  Total MCMC iterations per chain (default = 12500).

- n.burnin:

  Number of burn-in iterations (default = 20% of `n.iter`).

- dsep:

  Logical; if `TRUE`, performs d-separation tests to evaluate the global
  fit of the DAG against the data. Highly recommended for causal
  validation.

- family:

  A named character vector specifying the distribution for each response
  (e.g., `c(Y = "poisson", M = "gaussian")`). Supports "gaussian"
  (default), "poisson", "binomial", "gamma", "lognormal", "bernoulli",
  "ordinal", and "occupancy".

- latent:

  Optional character vector of latent (unmeasured) variables.

- parallel:

  Logical; if `TRUE`, runs MCMC chains in parallel.

- n.cores:

  Number of CPU cores for parallel execution.

- random:

  Optional formula for **global random intercepts** (e.g., `~(1|Site)`).
  Deprecated in favor of inline specification: `Y ~ X + (1|Site)`.

- levels:

  (Hierarchical) Optional named list mapping variables to their home
  levels. If omitted, `because` will attempt to auto-detect levels based
  on column availability.

- hierarchy:

  (Hierarchical) A string defining the nesting structure (e.g.,
  `"site > individual"`). Semicolons separate parallel branches (e.g.,
  `"site > obs; species > obs"`).

- link_vars:

  (Hierarchical) A named character vector specifying the columns used to
  link data frames across levels (e.g., `c(site = "SiteID")`).

- ...:

  Additional arguments passed to the underlying model engines.

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
# 1. Simple Path Model
df <- data.frame(Y = rnorm(100), M = rnorm(100), X = rnorm(100))
eqs <- list(M ~ X, Y ~ M + X)
fit <- because(eqs, data = df, dsep = TRUE)
summary(fit)

# 2. Hierarchical Model (Auto-detection)
# Suppose we have year-level data and individual-level data
year_df <- data.frame(Year = 1:5, Temp = rnorm(5))
ind_df  <- data.frame(Year = rep(1:5, each=10), Mass = rnorm(50))
data_list <- list(yr = year_df, ind = ind_df)

# Mass depends on Temp (cross-level link via Year)
fit_h <- because(
  equations = list(Mass ~ Temp),
  data      = data_list,
  hierarchy = "yr > ind",
  link_vars = c(yr = "Year")
)

# 3. Random Intercepts (lme4-style)
# Note: Grouping variables (Site) must be in the data
df$Site <- rep(1:10, each=10)
fit_re <- because(list(Y ~ X + (1|Site)), data = df)
} # }
```
