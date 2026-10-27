# Announcing because v1.8.0: Bayesian Estimation of Causal Effects in R

![because hex
sticker](https://raw.githubusercontent.com/because-pkg/because/master/man/figures/logo.png)

The **`because`** package is now publicly available in version
**1.8.0**, providing a unified R framework for specifying, fitting, and
evaluating Bayesian Structural Equation Models (BSEMs) with an explicit
focus on causal inference.

The package bridges Judea Pearl’s graphical causal framework (pearl,
2009) with applied Bayesian MCMC estimation, enabling researchers to
evaluate complex networks of observational hypotheses, test whether
causal graphs are supported by empirical data, and estimate direct and
indirect causal pathways across mixed distributions and structured
dependencies.

Observational studies in ecology, evolutionary biology, epidemiology,
and the social sciences routinely investigate systems with
interconnected causal pathways. While Structural Equation Modeling (SEM)
is a standard tool for such systems, conventional implementations often
face several challenges:

- **Non-Gaussian traits**: Real-world processes involve count data,
  zero-inflation, binary outcomes, and ordered categories that violate
  multivariate normal assumptions.
- **Complex dependencies**: Observational units often exhibit
  hierarchical clustering, phylogenetic relatedness, or spatial
  autocorrelation that cannot be addressed by naive linear models.
- **Missing data and measurement error**: Removing incomplete records
  introduces bias, while ignoring measurement uncertainty inflates
  false-positive rates.
- **Causal evaluation**: Traditional SEM is often interpreted as
  descriptive path analysis rather than an explicit causal framework
  with falsifiable conditional independence claims.

`because` addresses these challenges by translating standard R model
formulas into compiled Bayesian graphical models, providing formal tools
for causal model evaluation, diagnostic verification, and marginal
effect estimation.

------------------------------------------------------------------------

## Core Features

### 1. Multi-Engine Backend

Structural equations are specified using familiar R formula syntax.
`because` automatically compiles and executes the model using any of
three backends: \* **JAGS**: Robust Gibbs and slice sampling. \*
**NIMBLE**: High-performance compiled C++ MCMC with automated samplers
(including Elliptical Slice Sampling for structured covariance). \*
**NumPyro**: GPU- and CPU-accelerated Hamiltonian Monte Carlo (HMC/NUTS)
via Python’s JAX library, designed for large datasets.

### 2. Formal Causal Evaluation (d-Separation & m-Separation)

`because` tests whether the hypothesized causal graph is consistent with
the data. When `dsep = TRUE`, the package automatically identifies the
minimal basis set of conditional independencies implied by the graph
(Shipley, 2000) and performs fully Bayesian conditional independence
tests. To provide formal Bayesian evidence of equivalence rather than
relying on failure to reject a null hypothesis, conditional
independencies can be evaluated and visualized against a Region of
Practical Equivalence via `summary(..., rope = TRUE)` and
[`plot_rope()`](https://because-pkg.github.io/because/reference/plot_rope.md).
Maximal Ancestral Graphs (MAGs) are supported when unmeasured
confounding is present.

### 3. Average Marginal Effects on Response Scales (`marginal_effects()`)

In generalized structural equation models with non-Gaussian outcomes,
raw path coefficients are estimated on nonlinear link or latent scales
(such as log-odds for logistic paths or log-rates for Poisson paths),
making them difficult to interpret in terms of observable changes.
[`marginal_effects()`](https://because-pkg.github.io/because/reference/marginal_effects.md)
converts posterior parameter draws into Average Marginal Effects (AME)
on the **natural, observable scale of each respective response
variable** (e.g., probability differences $`\Delta P(Y = 1)`$ for binary
outcomes, or expected count differences $`\Delta \mathbb{E}[Y]`$ for
count responses). When paired with `plot_coef(type = "marginal")`, paths
are automatically grouped by their corresponding response variable,
enabling intuitive comparisons of effect sizes.

> **Standardized predictors in Generalized SEMs**: When continuous
> variables are z-standardized (mean = 0, SD = 1) prior to analysis,
> Gaussian-to-Gaussian path coefficients are already expressed in
> standard-deviation units ($`\beta^*`$). For non-Gaussian responses
> that cannot be standardized (e.g., integer counts or binary
> presence/absence),
> [`marginal_effects()`](https://because-pkg.github.io/because/reference/marginal_effects.md)
> translates a 1-SD shift in standardized predictors into tangible units
> of that response (such as expected count shifts
> $`\Delta \mathbb{E}[Y]`$ or absolute probability shifts
> $`\Delta P(Y = 1)`$), enabling direct comparison of predictor
> influence.

### 4. Automatic Missing Data Imputation & Measurement Error

Missing values in response and predictor variables are automatically
imputed under Missing at Random (MAR) or Missing Completely at Random
(MCAR) assumptions. Measurement error can be incorporated through
repeated measures or user-specified standard errors.

### 5. Broad Distribution Family Support

Supported likelihoods include: \* **Gaussian / Normal** (identity link)
\* **Binomial / Bernoulli** (logistic link) \* **Poisson** (log link) \*
**Negative Binomial** (overdispersed counts, log link) \*
**Zero-Inflated Poisson (ZIP)** and **Zero-Inflated Negative Binomial
(ZINB)** \* **Ordinal** (cumulative threshold model for ordered factors)
\* **Multinomial** (categorical choice models)

Additional custom families can be added via the S3 generic extension
system.

### 6. Integrated Bayesian Diagnostics Suite

`because` includes a dedicated diagnostic suite for evaluating model
assumptions: \* `check_model()`: Unified diagnostic visualizer. \*
`check_normality()` / `plot_qq()`: Randomized quantile residuals with
posterior predictive envelopes. \* `check_homoscedasticity()`: Binned
residual dispersion across fitted quantiles. \* `check_outliers()`:
Posterior predictive discrepancy detection. \* `check_collinearity()`:
Variance Inflation Factors (VIF) adapted for Bayesian models. \*
`check_residual_structure()`: Variance Partition Coefficients (VPC /
$`\lambda_{\text{Bayes}}`$) for phylogenetic or spatial random effects.
\*
[`pp_check()`](https://because-pkg.github.io/because/reference/pp_check.md):
Posterior predictive checks comparing observed data to simulated
replicated datasets.

------------------------------------------------------------------------

## Quick Example: Core `because` Workflow

The following self-contained example simulates a small ecological
dataset ($`n = 150`$), evaluates the causal graph via d-separation
(`dsep = TRUE`), fits the full simultaneous model (`dsep = FALSE`),
inspects the DAG with estimated posterior path coefficients, and
computes Average Marginal Effects:

``` r

library(because)

# 1. Simulate an ecological dataset
set.seed(42)
n <- 150
habitat   <- rnorm(n)
prey      <- rnorm(n, mean = 0.6 * habitat, sd = 0.7)
predators <- rnorm(n, mean = 0.5 * prey, sd = 0.7)
abundance <- rpois(n, lambda = exp(1.0 + 0.5 * prey - 0.4 * predators))
dat       <- data.frame(habitat, prey, predators, abundance)

# 2. Specify the causal model as a list of structural equations
eqs <- list(
  prey      ~ habitat,
  predators ~ prey,
  abundance ~ prey + predators
)

# 3. Step 1: Test the conditional independencies implied by the causal model using d-separation
dsep_res <- because(
  equations = eqs,
  data      = dat,
  family    = list(abundance = "poisson"),
  dsep      = TRUE
)
summary(dsep_res, rope = TRUE)
plot_rope(dsep_res)

# 4. Fit the full Bayesian SEM
fit_full <- because(
  equations = eqs,
  data      = dat,
  family    = list(abundance = "poisson"),
  dsep      = FALSE
)

# 5. Inspect parameter estimates and plot the fitted causal DAG
summary(fit_full)
plot_dag(fit_full)

# 6. Check model calibration via posterior predictive checks
pp_check(fit_full, resp = "abundance", type = "bars")

# 7. Compute Average Marginal Effects (AME) on each response's natural scale
ame <- marginal_effects(fit_full)
ame

# Or visualize effect sizes grouped by response variable
plot_coef(fit_full, type = "marginal")
```

------------------------------------------------------------------------

## Phylogenetic Comparative SEM: `because.phybase`

In comparative evolutionary biology species share phylogenetic history,
violating the assumption of independent and identically distributed
(i.i.d.) observations.

The **`because.phybase`** extension module integrates Phylogenetic
Bayesian Structural Equation Modeling (PhyBaSE) seamlessly into the
`because` workflow. Using an optimized random effects formulation
($`y_i = \mu_i + u_i + \epsilon_i`$, where
$`u \sim \mathcal{N}(0, \sigma_{\text{phylo}}^2 \mathbf{V})`$ and
$`\epsilon \sim \mathcal{N}(0, \sigma_{\text{res}}^2 \mathbf{I})`$),
Pagel’s $`\lambda`$ is estimated directly from variance components
($`\lambda = \sigma_{\text{phylo}}^2 / (\sigma_{\text{phylo}}^2 + \sigma_{\text{res}}^2)`$)
while pre-computing the inverse phylogenetic precision matrix once.

The example below demonstrates this using the classic Rhinogradentia
dataset and tree from the Phylogenetic Path Analysis chapter in *Modern
Phylogenetic Comparative Methods* (Gonzalez-Voyer & von Hardenberg,
2014) and von Hardenberg & Gonzalez-Voyer (2025):

``` r

library(because)
library(because.phybase)

# 1. Load the Rhinogradentia dataset and phylogenetic tree
data(rhino.dat, package = "because.phybase")
data(rhino.tree, package = "because.phybase")

rhino.dat[2:6] <-scale(rhino.dat[2:6])

# 2. Specify the causal model (Model 8 from Gonzalez-Voyer & von Hardenberg, 2014):
#    - Body mass (BM) affects litter size (LS)
#    - Body mass (BM) and range size (RS) affect nose length (NL)
#    - Nose length (NL) affects dispersal distance (DD)
sem8_eq <- list(
  LS ~ BM,
  NL ~ BM + RS,
  DD ~ NL
)

# 3. Test conditional independencies accounting for phylogenetic structure
dsep_phylo <- because(
  equations = sem8_eq,
  data      = rhino.dat,
  structure = rhino.tree,
  id_col    = "SP",
  dsep      = TRUE
)
summary(dsep_phylo)
plot_rope(dsep_phylo)

# 4. Fit the simultaneous phylogenetic SEM
fit_phylo <- because(
  equations = sem8_eq,
  data      = rhino.dat,
  structure = rhino.tree,
  id_col    = "SP",
  dsep      = FALSE
)

# 5. Summarize parameters (path coefficients, Pagel's lambda, variance components)
summary(fit_phylo)
plot_dag(fit_phylo)
```

`because.phybase` also supports: \* **Phylogenetic uncertainty**:
Integrating over posterior tree distributions (`multiPhylo` objects) to
propagate topological and branch-length uncertainty into causal effect
estimates. \* **Measurement error**: Handling repeated intraspecific
measurements across tips to separate evolutionary covariance from
within-species variance and sampling error.

------------------------------------------------------------------------

## The `because` Ecosystem

The package uses a modular S3 generic architecture that separates core
causal inference from domain-specific dependency structures:

- **`because` (Core)**: General Bayesian causal SEM, d-separation,
  missing data handling, marginal effects, and diagnostic suite.
- **`because.phybase`**: Phylogenetic comparative methods (PhyBaSE),
  incorporating single trees (`phylo`) and distributions of trees
  (`multiPhylo`).

------------------------------------------------------------------------

## Installation & Documentation

The packages can be installed from GitHub:

``` r

# Install core package
remotes::install_github("because-pkg/because@v1.8.0")

# Install phylogenetic extension
remotes::install_github("because-pkg/because.phybase")
```

- **Core Documentation**: <https://because-pkg.github.io/because/>
- **`because.phybase` Documentation**:
  <https://because-pkg.github.io/because.phybase/>
- **GitHub Repositories**:
  - <https://github.com/because-pkg/because>
  - <https://github.com/because-pkg/because.phybase>
- **Bug Reports & Issues**:
  <https://github.com/because-pkg/because/issues>

------------------------------------------------------------------------

## References

Gonzalez-Voyer, A., & von Hardenberg, A. (2014). An introduction to
phylogenetic path analysis. In L. Z. Garamszegi (Ed.), *Modern
Phylogenetic Comparative Methods and Their Application in Evolutionary
Biology: Concepts and Practice* (pp. 201–229). Springer, Berlin,
Heidelberg. <https://doi.org/10.1007/978-3-662-43550-2_8>

Pearl, J. (2009). *Causality: Models, Reasoning and Inference* (2nd
ed.). Cambridge University Press.
<https://doi.org/10.1017/CBO9780511803161>

Shipley, B. (2000). A new inferential test for path models based on
directed acyclic graphs. *Structural Equation Modeling*, 7(2), 206–218.
<https://doi.org/10.1207/S15328007SEM0702_4>

von Hardenberg, A., & Gonzalez-Voyer, A. (2025). PhyBaSE: A Bayesian
structural equation model approach to causal inference in phylogenetic
comparative analyses. *Methods in Ecology and Evolution*, 16(6),
1136–1148. <https://doi.org/10.1111/2041-210X.70044>
