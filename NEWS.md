# because 1.9.0 (2026-10-25)

## New Features & Enhancements
* **Posterior Predictive Summaries (`predict.because`)**: Added an S3 `predict()` method for `because` models returning observation-level posterior predictive summaries (`Estimate`, `Est.Error`, and credible quantiles) or raw posterior draws (`summary = FALSE`), complementing `fitted()` and `posterior_predict()`.
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

# because 1.8.0 (2026-09-24)

* **Bayesian Model Diagnostics Suite**: Added `check_model()`, `check_normality()`, `plot_qq()`, `check_homoscedasticity()`, `check_outliers()`, `check_collinearity()`, and `check_residual_structure()`.
* **Bayesian Residuals & Fitted Values**: Added `residuals.because()` (supporting `"quantile"` randomized PIT residuals, `"response"`, and `"pearson"`) and `fitted.because()`.
* **NumPyro Documentation & Installation**: Expanded documentation for the NumPyro backend and JAGS 4.x / 5.0 setup.

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
* **Counterfactual Interventions (`do()`)**: Added `because_do()` (`do()`) and `plot_counterfactual()` to simulate structural interventions (atomic, additive shifts, multipliers, and stochastic policies) using Pearl's do-calculus.
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
