# Changelog

All notable changes to `because` are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/2.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

# because 1.9.0 (2026-11-11)

## New Features & Improvements
* **Standardization on Latent Liability and Common SD Scale (`standardize()`)**:
  * Added S3 generic and method `standardize()` to rescale path coefficients from fitted models and d-separation test objects onto a common, dimensionless standard deviation metric.
  * **Continuous Gaussian Traits**: Scaled by the ratio of standard deviations $\beta^* = \beta \times (\text{SD}(X) / \text{SD}(Y))$ (preserving scale if variables were standardized prior to fitting).
  * **Binary / Categorical Responses (Binomial/Bernoulli/Ordinal Logit)**: Raw log-odds coefficients are transformed *a posteriori* onto the continuous latent liability scale following [McKelvey & Zavoina (1975)](https://doi.org/10.1080/0022250X.1975.9989847), [Menard (2004)](https://doi.org/10.1198/000313004X946), and [Nakagawa et al. (2017)](https://doi.org/10.1098/rsif.2017.0213):
    $$\beta^* = \frac{\beta \times \text{SD}(X)}{\sqrt{\sum \sigma_{\text{random}}^2 + \pi^2/3}}$$
    where $\pi^2/3 \approx 3.289868$ is the theoretical disturbance variance of the standard logistic link. Full posterior uncertainty of random effects is propagated draw-by-draw into $\beta^*$.
  * **Count Responses (Poisson / Negative Binomial)**: Added full standardization support for count responses on the log rate scale ($\ln \mu$) across three methods:
    * `method = "error"` (default): $\beta^* = \frac{\beta \times \text{SD}(X)}{\sqrt{\sum \sigma_{\text{random}}^2 + \ln(1 + 1/\bar{Y})}}$ following the GLMM latent distribution variance framework of [Nakagawa et al. (2017)](https://doi.org/10.1098/rsif.2017.0213).
    * `method = "total"`: Full total variance incorporating the empirical variance of the fixed linear predictor $\text{Var}(\hat{\eta}_{\text{fixed}})$ alongside random effects and disturbance variance.
    * `method = "mplus"`: Standardizes by the standard deviation of the linear predictions $s_{xb} = \sqrt{\text{Var}(\hat{\eta}_{\text{fixed}}) + \sum \sigma_{\text{random}}^2}$, exactly reproducing Mplus `StdYX` for Poisson regression ([UCLA OARC FAQ](https://stats.oarc.ucla.edu/mplus/faq/how-does-mplus-calculate-the-standardized-coefficients-based-on-a-poisson-model/)).
  * **Semi-Standardization for Binary Predictors (`binary_predictors = "raw"`, default)**: Following [Gelman (2008)](https://doi.org/10.1002/sim.3107), [Bollen (1989)](https://doi.org/10.1002/9781118619179), and [Grace & Bollen (2005)](https://doi.org/10.1890/04-0659), binary dummy predictors ($X \in \{0, 1\}$) are kept on their natural $0 \to 1$ step while standardizing the outcome, so that $\beta^*$ directly represents **Cohen's $d$** (group difference in standard deviations of the outcome). Optional `binary_predictors = "sd"` supports classical full standardization ($s_X = \sqrt{p(1-p)}$; [Menard 2004](https://doi.org/10.1198/000313004X946)).
  * **Automated Binary Detection Fallback**: Automatically infers binomial links from $\{0, 1\}$ data even if the `family` attribute was omitted from the fitted object.
  * **High-Efficiency In-Place Updates**: Selectively updates only standardized parameter rows in `summary$statistics` and `summary$quantiles`, preserving precomputed MCMC convergence diagnostics (`Rhat`, `Rhat_max`, `n.eff`) and avoiding memory duplication of large ensemble posterior arrays.
* **ROPE and Ridge Visualization Harmonization (`plot_rope()`, `summary()`)**:
  * When evaluating standardized models (`standardized = TRUE`), `resolve_rope_for_param()` automatically applies the uniform standard deviation benchmark `c(-0.10, 0.10)` across all paths and d-separation tests.
  * `plot_rope()` plots standardized regression slopes ($\beta^*$) and correlation coefficients ($\rho \in [-1, 1]$) together on a unified scale against the $\pm 0.10$ ROPE band.
  * `summary()` appends an informative footnote documenting standardized and latent liability evaluation.
* **Vectorized Diagnostics & Fast Ensemble Summaries**:
  * Vectorized `coda::effectiveSize()` and `coda::gelman.diag()` across candidate trees in `summary.because()`, restricting recalculation strictly to structural parameters (`alpha`, `beta`, `sigma`, `rho`, `lambda`).
  * Summaries of large multiPhylo ensembles (30+ trees, 60+ chains) now execute in under 0.1 seconds rather than minutes.
  * Primary summary tables now present within-tree $\hat{R}$ and worst-case tree $\hat{R}_{\text{max}}$, leaving between-tree sensitivity ($\hat{R}_{\text{tree}}$) to `check_multicov()`.
* **Memory-Safe WAIC & *A Posteriori* Evaluation (`because_waic()`)**:
  * `because_waic()` now returns precomputed `model$WAIC` directly when available, and supports manual *a posteriori* evaluation on fitted NumPyro ensemble models without requiring MCMC refitting.
  * Wrapped WAIC formatting in exception handling in `run_numpyro_pipeline()` to ensure formatting or memory issues issue a warning rather than discarding completed MCMC samples.
* **Distributional Parsimony Rule for d-Separation**:
  * When testing auxiliary conditional independences between continuous and discrete/categorical variables, the continuous variable is prioritized as the regression response. This routes the test to fast ANOVA / linear mixed models, bypassing slow discrete links while preserving exact d-separation topology.
* **Liability-Scale Pagel's Lambda for Binary Traits**:
  * In JAGS models with phylogenetic structures and binomial traits, Pagel's $\lambda$ is now evaluated on the continuous latent liability scale using the theoretical $\pi^2/3$ disturbance variance.
* **Universal Multi-Structure Covariance Ensembles**: `structure` now natively accepts direct named lists of candidate covariance structures across any domain (e.g., phylogenetic backbones `list(Jetz = jetz_trees, McTavish = mctav_trees)`, spatial dispersal kernels, genomic relationship matrices, or animal social networks). When a named list of structures or collections is supplied, `because` automatically flattens candidate structures into a unified collection, assigns group/backbone attributes, and partitions chains across candidate structures without requiring explicit control wrappers.
* **Universal Multi-Structure Convergence & Sensitivity Diagnostic (`check_multicov()`)**: Added `check_multicov()` to evaluate MCMC convergence and structural sensitivity for any model fitted with multiple candidate covariance structures across all engines (NumPyro, JAGS, and NIMBLE). Computes a 3-tier variance decomposition:
  * $\hat{R}_{\text{MCMC}}$: Within-structure MCMC convergence across replicate chains ($\max_k \hat{R}_k$).
  * $\hat{R}_{\text{within\_group}}$: Parameter sensitivity across candidate structures within the same group, hypothesis, or backbone pseudo-posterior.
  * $\hat{R}_{\text{between\_group}}$: Sensitivity and divergence between competing groups, kernels, or backbone hypotheses using Gelman-Rubin ANOVA decomposition.
  * **Exact Law of Total Variance Uncertainty Decomposition**: Orthogonally partitions total posterior uncertainty for each parameter into additive components that sum to exactly 100%: `% MCMC` (data sampling/estimation uncertainty conditional on a structure), `% Structure` (branch-length / dating / local topological variation within group), and `% Group` (deep structural conflict between competing hypotheses/backbones).
  * Informative diagnosis classification: `"Invariant"` (robust across all structures and groups), `"Structure-Sensitive (within group)"`, `"Group-Sensitive"`, or `"Check MCMC"`.
  * Visual diagnostic plots: 2-panel publication figure displaying the $\hat{R}$ convergence/sensitivity forest plot (left) alongside a 100% stacked horizontal bar plot of the uncertainty variance decomposition (right), aligned row-for-row per parameter. Trace and frequency plots are generated for discrete jump sampling (`method = "dcat"`).
* **Structure Sensitivity Reporting in d-Separation (`Rhat_struct`)**: Both `dsep_numpyro.R` and `dsep_dispatch.R` now report structural sensitivity alongside within-structure $\hat{R}$ across auxiliary conditional independence regressions.
* **Unified `n.chains` Parameter Across Multi-Structure Engines**: Standardized `n.chains` across all structure handlers, model compiler routines, and documentation.
* **Preserved Structure and Group Metadata on Fitted Models**: Fitted `because` model objects in NumPyro, JAGS, and NIMBLE now consistently preserve `fit$structure`, `fit$structures`, and `fit$input$structure` with their attached metadata attributes (`"group"`, `"backbone"`, `"tree_method"`, `"n.chains"`).
* **Multi-Response Q-Q Envelope Diagnostics (`plot_qq()`)**: Calling `plot_qq(fit)` or `check_normality(fit)` on a fitted SEM now automatically checks all endogenous response variables when `resp = NULL`, rendering a multi-panel faceted grid with empirical 50% (interquartile) and 95% credible simulation envelopes ([Landwehr et al. 1984](https://doi.org/10.1080/01621459.1984.10478048)). Evaluates randomized quantile residuals via the probability integral transform (PIT against $y^{\text{rep}}$, [Dunn & Smyth 1996](https://doi.org/10.1080/10618600.1996.10474708)) with posterior predictive skewness and kurtosis realized discrepancies ([Gelman et al. 1996](https://www.jstor.org/stable/24306036)). Single responses (e.g. `resp = "Y1"`) or custom subsets (e.g. `resp = c("Y1", "Y2")`) remain fully supported.
* **Conditional Posterior Predictive Checks for SEM Latent Measurement Blocks (`conditional = "auto"`)**: In structural equation models with unmeasured latent variables or induced residual covariances (e.g., `Y1 ~ Factor`, `Y2 ~ Factor` via the MAG framework), `posterior_predict()` and `pp_check()` now automatically condition on correlated companion traits using the closed-form conditional multivariate Gaussian distribution ($y_1 \mid \mathbf{y}_C^{\text{obs}}$). This resolves the issue where marginal root-node simulations produced an artificial single bell curve against bimodal or multi-modal trait distributions. Controlled via `conditional = "auto"` (default), `conditional = TRUE`, or `conditional = FALSE` (pure unconditional forward causal generation), with optional explicit `condition_on`.
* **Informative Conditioning Labels in `pp_check()`**: When conditional predictions are active, `pp_check()` automatically documents the conditioning traits in the plot subtitle (e.g. `Type: dens_overlay | draws: 50 | Conditional on: Y2, Y3`), making diagnostic displays transparent for publication and reviewer evaluation.
* **KaTeX Math Rendering in Reference Documentation**: Added `math-rendering: katex` to `_pkgdown.yml`, ensuring LaTeX equations in roxygen references (such as mathematical formulas for quantile residuals, VIF thresholds, and variance decompositions) render cleanly in online documentation.

## Bug Fixes
* **Exact Attribute Matching for Structure Dimensions**: Fixed a partial matching bug in `process_because_structures()` where base R `attr(structure_obj, "n")` inadvertently matched `attr(structure_obj, "n.chains")` because `exact = TRUE` was omitted, causing dimension mismatch errors when fitting multi-chain ensembles.
* **Safe Unpacking of Multi-Array Structure Data**: In `process_because_structures()`, structure data returned as nested lists of matrices/arrays (e.g., `list(eigvals = ..., eigvecs = ...)` in NumPyro) is now unpacked when inferring structure observation dimensions $N$.
* **Preserved Variance for Marginalized Random Effects in Posterior Predictions**: When species-level phylogenetic or hierarchical random effects are analytically marginalized during MCMC (such as `dmnorm` or spectral marginalization) without storing individual species vectors $u_i$, `posterior_predict()` now automatically incorporates their sampled variance ($\sigma_{\text{phylo}}^2$) into predictive simulation noise instead of discarding it. This prevents `pp_check()` and `plot_qq()` from collapsing into artificial zero-width needles at $x = 0$.
* **Roxygen Topic Association for Residuals**: Moved internal helper `%||%` before the roxygen comment block in `R/residuals.R`, ensuring that the reference index correctly documents `residuals.because` and `fitted.because` rather than exposing the internal operator.

---

# because 1.8.3 (2026-11-03)

## New Features & Improvements
* **NumPyro ROPE and Ridge Plots for d-Separation**: NumPyro d-separation runs now preserve raw posterior draws across individual conditional independence tests in `dsep_results`. This enables full ROPE support for the NumPyro engine, including `P_ROPE` columns in `summary(fit, rope = ...)` and graphical ridge visualizations with `plot_rope(fit)`.
* **Sensible Weakly Informative Prior for Negative Binomial Dispersion (`r`)**: Replaced the legacy BUGS default prior `r ~ dgamma(0.01, 0.01)` (which concentrated 96% of its prior mass below 1 and had a median of $10^{-29}$, unintentionally placing extreme weight on astronomical overdispersion) with a weakly informative `dgamma(2, 0.1)` prior (mean = 20, variance = 200, aligning with modern recommendations and the NumPyro specification in `because_py`).
* **Custom Dispersion and Zero-Inflation Priors**: The model compiler now routes Negative Binomial dispersion parameters (`r_<var>`) and zero-inflation probabilities (`psi_<var>`) through `get_prior()`, allowing users to provide custom priors either by parameter name (e.g. `priors = list(r_Y = "dgamma(5, 0.5)")`, `psi_Y = "dbeta(2, 5)"`) or by generic parameter type (`r = "dgamma(5, 0.5)"`, `psi = "dbeta(1, 1)"`).
* **Automated Prior Scale Diagnostic Check**: Added `check_prior_scale()` which inspects the empirical variance of continuous Gaussian variables before model compilation. If unstandardized variables with large standard deviations ($\text{SD} > 20$) are detected without user-specified custom priors, `because()` issues an informative warning alerting the user that default priors (e.g. `sigma ~ dunif(0, 100)`, $\alpha, \beta \sim \text{Normal}(0, 100)$) assume standardized variables ($\text{mean} = 0, \text{SD} = 1$), advising them to standardize continuous traits or supply scale-appropriate custom priors to prevent unintended prior truncation. The check is non-blocking (the model runs anyway) and strictly restricted to continuous Gaussian responses (ignoring Binomial, Poisson, and other discrete/count families).
* **Direct Custom Prior Overrides on Residual Standard Deviation (`sigma_*_res`, `sigma_e_*`)**: `get_precision_prior()` now natively recognizes custom priors specified directly on the standard deviation scale for Gaussian residuals (e.g. `sigma_Y_res = "dunif(0, 500)"` or `sigma_e_Y = "dunif(0, 500)"`), automatically translating them into the corresponding JAGS precision statements (`tau_res_Y <- 1 / (sigma_Y_res * sigma_Y_res)`).
* **Documentation & Vignette Clarifications on Prior Scaling**: Updated `vignettes/07_custom_priors.Rmd` with dedicated guidance on variable scaling and prior bounds, clarifying that default priors are designed for standardized traits to ensure comparable effect sizes, and illustrating how to specify custom priors for unstandardized raw data.
* **Expanded Test Suite**: Added `tests/testthat/test-check_prior_scale.R` and `tests/testthat/test-dispersion-priors.R` covering scale diagnostics, custom prior overrides, discrete family exemptions, quiet mode, and dispersion/zero-inflation custom priors, bringing the package test suite to 35 test files and 296 passing tests.

## Bug Fixes
* **Fixed `because_continue()` Parallel Rebuild Data Argument**: Fixed an issue in `because_continue()` where Strategy 2 (parallel model rebuild across worker nodes) referenced an unbound variable `data` instead of `object$data`. When worker nodes evaluated the closure, `data` resolved to the base R function `utils::data`, causing JAGS compilation to crash with `Error in checkForRemoteErrors(val) : 3 nodes produced errors; first error: data must be a list or environment`. `because_continue()` now properly defines and passes `model_data <- object$data`, safely checks whether the live JAGS external pointer is valid via `object$model$iter()` before attempting in-session continuation, and falls back cleanly to parallel or sequential warm-start rebuild when reading saved `.rds` / `.Rdat` model files or continuing across separate R sessions. Added comprehensive regression tests in `test-because_continue.R`.
* **Restored Observation Precision & Conjugacy for Latent MAG Indicators in JAGS**:
  * Fixed an issue in `generate_likelihoods()` where measurement indicator variables with correlated residual errors (MAG approach, e.g. `Brain` and `Mass` for latent `Size`) were compiled using `tau_res_<var>` instead of `tau_obs_<var>`. This had caused the model compiler to assign an unidentifiable non-conjugate `dunif(0, 100)` prior on residual standard deviation (`sigma_<var>_res`), creating four competing variance components (`TAU_res_*` Wishart matrix, `sigma_<var>_res`, `sigma_obs_*`, and `sigma_<struct>_*`) that led to boundary collapse (`sigma_res -> 0`), multimodal variance swapping, and severe MCMC mixing / convergence issues (elevated $\hat{R} > 2.0$ and low ESS).
  * Restored conjugate observation precision `tau_obs_<var> ~ dgamma(1, 1)` and derived standard deviation `sigma_obs_<var> <- 1/sqrt(tau_obs_<var>)`, eliminating the conflicting `sigma_<var>_res` prior.
  * Eliminated duplicate phantom prior declarations (`tau_u_<var>_<struct>` vs `tau_u_<struct>_<var>`) in `generate_model_priors()`.
  * Corrected parameter regex filtering in `summary.because()` and `prepare_monitors_and_inits()` so redundant structural and residual variances for MAG indicators are omitted from monitoring and summary tables while retaining observation and residual error parameters (`sigma_obs_*`, `sigma_res_*`, and `rho_*`).
* **Multiscale Hierarchical d-Separation & Cross-Level Indexing in NumPyro**: Fixed three issues in the NumPyro backend (`because` and `because_py`) for multiscale / hierarchical models:
  * `run_dsep_numpyro()` now injects applicable structural random effects (`+ (1 | <structure>)`) via `compile_equation_string` into each d-separation sub-model and initializes structure index vectors using `hierarchical_info$structure_levels`.
  * `because` now passes explicit level-to-variable and structure mappings (`levels`) to `because_py`, preventing cross-level index collisions (`<level>_idx_<child>`) when two orthogonal hierarchy branches have identical sample sizes (e.g. `N_site == N_species`) across fixed predictors, `I()` deterministic terms, and random effects.
  * `prepare_structure_data.matrix()` now preserves the covariance matrix (`VCV`) alongside the precision matrix (`Prec`) so NumPyro computes the Cholesky factor of the covariance matrix rather than its inverse.
* **Partial Named Vector `family` Handling in d-Separation and `plot_rope()`**: Fixed a `"subscript out of bounds"` error when `family` is passed as a partial named character vector (e.g. `family = c(Abundance = "poisson")`) during d-separation testing, and updated `plot_rope()` to apply response-specific family ROPE thresholds (`[-0.05, 0.05]` for Poisson/NegBin, `[-0.18, 0.18]` for Binomial/Ordinal, `[-0.10, 0.10]` for Gaussian) per test in multi-family models.
* **Structured Random Effects for Non-Gaussian Responses**: Fixed a regression (introduced when the model compiler was modularized) where responses with a link function (e.g. `poisson`, `binomial`, `negbinomial`) combined with a `structure` (phylogenetic, spatial or group) failed JAGS compilation with `Unknown variable sigma_<structure>_<response>`. The structural variance priors (`tau_u_*`, `sigma_*`) are now declared for all families, while residual variance (`tau_res_*`) and Pagel's $\lambda$ remain restricted to families with additive residual variance. Added a regression test in `test-distributions.R`.
* **Latent MAG Indicators Excluded from Pagel's Lambda and Redundant Variance Nodes**: Measurement indicator variables linked to marginalized latent variables in the Mixed Ancestral Graph (MAG) approach (e.g. `Brain`, `HWI`, `Mass` for latent `Size`) have induced residual covariance structures (`rho_*` / `TAU_res_*`), for which a univariate Pagel's $\lambda$ is undefined. These indicator variables are now systematically excluded from Pagel's $\lambda$ calculations and monitoring across both `because` and `because_py` (v0.2.10). In addition, redundant internal variance partitioning nodes (`sigma_tot_*` and `sigma_*_<structure>`) for these indicators are excluded from model reporting, matching the clean output of the JAGS engine.
* **Robust Monitor Regex Matching for Structural Group Suffixes**: Strengthened regex patterns in `prepare_monitors_and_inits()` so parameters with structural suffixes (such as `lambda_<var>_<structure>` and `alpha_<var>_<structure>`) are correctly matched against response variables and filtered out for MAG exogenous variables.

---

# because 1.8.2 (2026-10-29)

## Bug Fixes & Improvements
* **NumPyro Latent MAG Induced Correlation Reporting**: Fixed an issue where induced correlation parameters (`rho_*`) resulting from marginalizing latent variables (MAG approach) were missing from `summary()` when using `engine = "numpyro"`. `because` now forwards `induced_cors` into `run_numpyro_pipeline()` and `because_py` (v0.2.9), and extracts the off-diagonal correlation coefficient directly from the 4D LKJ Cholesky factor matrix (`L_corr_*[..., 2, 1]`) to register it as `rho_*`.
* **Hierarchical Random Effects Initialization**: Fixed an uninitialized variable error (`object 'row_ids' not found`) in `prepare_random_effects_data()` when fitting multi-level / hierarchical models.
* **Plot Label Formatting**: Normalized whitespace handling in `plot_dsep()` to ensure clean formatting for tests with empty conditioning sets.

---

# because 1.8.1 (2026-10-28)

## Bug Fixes & Improvements
* **Categorical Dummy Variable Name Sanitization**: Fixed JAGS/NIMBLE compilation errors when categorical variables or their levels contain periods (`.`) by systematically sanitizing base column and level names (e.g. `Sample.size` $\rightarrow$ `Sample_size_medium`).
* **Categorical Dummy Disambiguation in Plots**: `plot_rope()` and `plot_dsep()` now disambiguate individual dummy coefficients for categorical predictors (e.g. `Brain _||_ Sample.size (medium)`), preventing `plot_rope()` from unintentionally pooling all category draws into a single merged density.
* **Observation Row Ordering in Structure Preparation**: Retained observation `row_ids` in flat models during random effects / structure data formatting to prevent phylogenetic tip alignment issues.

---

# because 1.8.0 (2026-10-26)

## New Features & Enhancements
* **Multi-Family & Principled ROPE Evaluation**: Added distribution-aware Region of Practical Equivalence (ROPE) evaluation to `summary()` and `plot_rope()`. Supports customized lists by family or response, uses principled defaults (Gaussian: `[-0.10, 0.10]`, Poisson/NegBin: `[-0.05, 0.05]`, Binomial/Ordinal: `[-0.18, 0.18]`), and displays human-readable rate/odds ratio interpretations alongside warnings for unstandardized continuous variables.
* **Bayesian Model Diagnostics Suite**: Added a comprehensive suite of model diagnostic tools:
  * `check_model()`: Unified diagnostic visualizer.
  * `check_normality()` & `plot_qq()`: Randomized quantile residual Q-Q plots with empirical predictive envelopes.
  * `check_homoscedasticity()`: Binned posterior residual standard deviations across fitted value quantiles.
  * `check_outliers()`: Posterior predictive discrepancy detection for extreme values.
  * `check_collinearity()`: Variance Inflation Factors (VIF) for Bayesian predictor structures.
  * `check_residual_structure()`: Posterior Variance Partition Coefficients (VPC / $\lambda_{\text{Bayes}}$) for phylogenetic, spatial, and custom grouping random effects.
* **Bayesian Residuals & Prediction Methods**:
  * Added `residuals.because()` supporting `"quantile"` (randomized PIT residuals; [Dunn & Smyth 1996](https://doi.org/10.1080/10618600.1996.10474708)), `"response"`, and `"pearson"` residuals.
  * Added `fitted.because()` to compute posterior expected values across Gaussian, Binomial, Poisson, Negative Binomial, ZIP, ZINB, Gamma, and Lognormal models.
  * Added `predict.because()` returning observation-level posterior predictive summaries (`Estimate`, `Est.Error`, and credible quantiles) or full posterior prediction arrays (`summary = FALSE`).
* **Expanded Posterior Predictive Checks (`pp_check`)**:
  * Added support for discrete and grouped `bayesplot` check types: `"bars"`, `"bars_grouped"`, `"stat_2d"`, `"intervals"`, `"ecdf_overlay"`, and `"scatter_avg"`, alongside `"dens_overlay"`, `"hist"`, `"boxplot"`, and `"stat"`.
  * Added the `group` argument for grouped posterior predictive checks.
  * Changed default `trim = FALSE` so extreme tail replications remain visible by default.
* **Generative Simulation for Additional Families**: `posterior_predict()`, `fitted()`, and `marginal_effects()` now support Zero-Inflated Poisson (`zip`), Zero-Inflated Negative Binomial (`zinb`), and `ordinal` families, as well as interaction terms (`:` / `_x_`) and cross-scale hierarchical predictor expansion.
* **Fully Bayesian Cross-Scale d-Separation**: Replaced frequentist PGLS cross-scale d-separation shortcuts with full Bayesian hierarchical aggregation across scales.
* **Engine-Aware MCMC Defaults**: `because()` now sets engine-appropriate default iteration, burn-in, and adaptation counts (`n.iter`, `n.burnin`, `n.adapt`) for `"jags"`, `"nimble"`, and `"numpyro"`.
* **Generalized Extension Architecture**:
  * Extended `dag_expand_hook()` to support `compound_groups` metadata for multi-parameter observational families.
  * Added `dsep_structure_hook()` to delegate structure subsetting during d-separation tests cleanly to extension packages (`because.phybase`).
  * Generalized NIMBLE's Elliptical Slice Sampler (ESS) assignment to cover spatial random effects alongside phylogenetic structures.

## Bug Fixes & Statistical Rigor
* **Synchronized Missing Data (`NA`) Alignment in Diagnostics**: Fixed an index-alignment bug in `pp_check()`, `residuals()`, and `check_normality()` where missing values (`NA`) in observed responses were stripped without subsetting the corresponding columns of `yrep`.
* **GLM Pearson Residuals**: `residuals(..., type = "pearson")` now standardizes by the family variance function $\sqrt{V(\mu)}$ for Poisson, Binomial, Negative Binomial, ZIP, and ZINB models instead of assuming Gaussian residual standard deviation.
* **Negative Binomial Dispersion in `posterior_predict()`**: Fixed parameter lookup to extract `r_<resp>` from posterior samples rather than defaulting to $r = 1$.
* **Response-Scale Marginal Effects for Count Models**: `marginal_effects()` now applies the exponential/zero-inflated inverse link for `poisson`, `negbinomial`, `gamma`, `lognormal`, `zip`, and `zinb` responses.
* **Random Number Generator (RNG) Hygiene**: Removed internal `set.seed()` calls from post-estimation functions (`posterior_predict()`, `fitted()`, `residuals()`, `plot_coef()`, `dsep`) so user `.Random.seed` state is never overwritten, and derived parallel worker chain seeds dynamically from R's RNG state.

---

# because 1.7.0 (2026-09-21)

* **Modular Multi-Engine Architecture**: Refactored the JAGS, NIMBLE, and NumPyro execution pipelines into dedicated backend modules (`backend_jags.R`, `backend_nimble.R`, `backend_numpyro.R`) and split the model compiler, d-separation engine, hierarchical data processor, and DAG plotting utilities into focused modules.
* **Incremental Model Refitting**: Added `update.because()` and `because_refit()` to support adding or removing equations while reusing cached d-separation tests.
* **Centered Random Effects**: Improved MCMC mixing for hierarchical random intercepts and made engine progress messages dynamic across JAGS, NIMBLE, and NumPyro.

---

# because 1.6.0 (2026-07-03)

* **NumPyro NUTS Controls**: Added `max_treedepth` and `adapt_delta` tuning controls for NumPyro's No-U-Turn Sampler (NUTS) and updated compatibility to `because_py` v0.2.8.
* **Prior & Hierarchy Refinements**: Tightened default slope priors for non-identity link functions to prevent numerical overflow during warmup, and improved multiscale hierarchy validation and level merging.
* **Exogenous Structural Priors**: Fixed structural prior generation for exogenous correlated variables in MAG models.

---

# because 1.5.0 (2026-06-25)

* **Model Call Preservation & Documentation**: Preserved the original matched call in fitted `because` objects to support incremental model updating (`update()`) and standardized roxygen documentation and cross-links across post-estimation and plotting functions.

---

# because 1.4.0 (2026-06-24)

* **Region of Practical Equivalence (`plot_rope`)**: Added `plot_rope()` and `rope` support in `summary.because()` for Bayesian equivalence testing, alongside customizable credible interval probabilities (`prob`).
* **MAG Orientation & m-Separation**: Improved Maximal Ancestral Graph (MAG) orientation rules, automatic latent variable detection, and m-separation basis set construction when unmeasured confounders (`induced_correlations`) are present.
* **Automatic Reticulate Environment Binding**: `.onLoad()` now automatically binds `reticulate` to the `because_env` virtual/conda environment before Python initialization.

---

# because 1.3.0 (2026-06-09)

* **NumPyro Engine Integration**: Added `engine = "numpyro"` and `install_because_numpyro()` for GPU/CPU-accelerated Hamiltonian Monte Carlo (NUTS) sampling via JAX and `because_py`.
* **Continuation Sampling (`because_continue`)**: Added `because_continue()` to extend MCMC chains for fitted JAGS and NIMBLE models without discarding existing burn-in and posterior samples.
* **Counterfactual Interventions (`do()`)**: Added `because_do()` (`do()`) and `summary.because_counterfactual()` to simulate structural interventions (atomic, additive shifts, multipliers, and stochastic policies) using Pearl's do-calculus.
* **Average Marginal Effects (`marginal_effects`)**: Added `marginal_effects()`, `plot_coef()`, and marginal-effect edge scaling in `plot_dag()` to compare paths across Gaussian, Binomial, Ordinal, and Multinomial responses on a unified expected-change scale.
* **Multiscale & Hierarchical Modeling**: Renamed `hierarchical` to `multiscale` with support for nested and cross-classified `(1 | group)` random effects and automatic scale inference.

---

# because 1.2.0 (2026-03-24)

* **Expanded Distribution Families**: Added Zero-Inflated Poisson (`zip`), Zero-Inflated Negative Binomial (`zinb`), and Negative Binomial (`negbinomial`) response families with NIMBLE custom distribution optimizations.
* **S3 Extension Architecture**: Introduced S3 generic hooks (`jags_structure_definition`, `prepare_structure_data`, `jags_family_likelihood`, `nimble_family_optimization`) allowing external packages (such as `because.phybase`) to register custom covariance structures and likelihood families.
* **d-Separation Visualization (`plot_dsep`)**: Added `plot_dsep()` for visualizing posterior conditional independence test estimates and credible intervals.

---

# because 1.1.0 (2026-03-14)

* **NIMBLE Engine Support**: Added `engine = "nimble"` for compiled C++ MCMC sampling, including automated Elliptical Slice Sampling (ESS) and block RW samplers for correlated structural effects.
* **Multilevel & Latent Variable Foundations**: Added hierarchical loop bound resolution, precision priors, and initial Maximal Ancestral Graph (MAG) induced correlation support.

---

# because 1.0.0 (2026-02-22)

* **First Stable Release of `because`**:
  * Automated Bayesian SEM specification and BUGS model compilation from lists of R formulas.
  * Native Bayesian d-separation testing (`dsep = TRUE`), mediation analysis (`because_mediation()`), WAIC (`because_waic()`), LOO (`because_loo()`), and multi-model comparison (`because_compare()`).
  * Support for deterministic nodes (`I()`), polynomial terms, categorical predictors, measurement error (`se` and `reps`), and automatic MAR/MCAR missing data imputation.
