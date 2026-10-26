# Package index

## Model Fitting & Estimation

Core functions for specifying and fitting Bayesian Structural Equation
Models.

- [`because()`](https://because-pkg.github.io/because/reference/because.md)
  : Run a Bayesian Structural Equation Model (Because)
- [`because_model()`](https://because-pkg.github.io/because/reference/because_model.md)
  : Generate BUGS Model Code and Parameter Map for Bayesian SEM
  (Because)
- [`because_continue()`](https://because-pkg.github.io/because/reference/because_continue.md)
  : Continue MCMC Sampling for a Fitted Model
- [`update(`*`<because>`*`)`](https://because-pkg.github.io/because/reference/update.because.md)
  : Update a Fitted because Model

## Causal Evaluation & d-Separation

Tools for evaluating causal DAG consistency and testing conditional
independencies.

- [`because_dsep()`](https://because-pkg.github.io/because/reference/because_dsep.md)
  : Extract d-separation statements from a structural equation model
- [`plot_dsep()`](https://because-pkg.github.io/because/reference/plot_dsep.md)
  : plot_dsep
- [`plot_rope()`](https://because-pkg.github.io/because/reference/plot_rope.md)
  : Plot Region of Practical Equivalence (ROPE)

## Causal Effects & Interventions

Functions for estimating marginal effects, causal mediation, and
structural interventions.

- [`marginal_effects()`](https://because-pkg.github.io/because/reference/marginal_effects.md)
  : Calculate Marginal Effects for Because Fit
- [`because_mediation()`](https://because-pkg.github.io/because/reference/because_mediation.md)
  : Automated Mediation Analysis
- [`do()`](https://because-pkg.github.io/because/reference/do.md)
  [`because_do()`](https://because-pkg.github.io/because/reference/do.md)
  [`summary(`*`<because_counterfactual>`*`)`](https://because-pkg.github.io/because/reference/do.md)
  : Causal Interventions for Structural Equation Models

## Model Diagnostics & Comparison

Posterior predictive checks, information criteria, and prediction
methods.

- [`summary(`*`<because>`*`)`](https://because-pkg.github.io/because/reference/summary.because.md)
  : Summary for Because Model
- [`pp_check()`](https://because-pkg.github.io/because/reference/pp_check.md)
  : Posterior Predictive Checks
- [`predict(`*`<because>`*`)`](https://because-pkg.github.io/because/reference/predict.because.md)
  : Predict Method for because Model Fits
- [`posterior_predict()`](https://because-pkg.github.io/because/reference/posterior_predict.md)
  : Posterior Predictive Samples
- [`extract_imputed()`](https://because-pkg.github.io/because/reference/extract_imputed.md)
  : Extract Imputed Values for Missing Data
- [`because_compare()`](https://because-pkg.github.io/because/reference/because_compare.md)
  : Compare Because Models
- [`because_loo()`](https://because-pkg.github.io/because/reference/because_loo.md)
  : Calculate LOO-CV for a Because Model
- [`because_loo_compare()`](https://because-pkg.github.io/because/reference/because_loo_compare.md)
  : Compare Models Using LOO-CV
- [`because_waic()`](https://because-pkg.github.io/because/reference/because_waic.md)
  : Calculate WAIC with Standard Errors for a Because Model

## Visualization

Plotting methods for DAGs, coefficients, and posterior distributions.

- [`plot_dag()`](https://because-pkg.github.io/because/reference/plot_dag.md)
  : Plot DAG from Equations or Fitted Model
- [`plot_coef()`](https://because-pkg.github.io/because/reference/plot_coef.md)
  : Plot Path Coefficients
- [`plot_posterior()`](https://because-pkg.github.io/because/reference/plot_posterior.md)
  : Plot Posterior Distributions for Model Comparison
- [`plot_dsep()`](https://because-pkg.github.io/because/reference/plot_dsep.md)
  : plot_dsep
- [`plot_rope()`](https://because-pkg.github.io/because/reference/plot_rope.md)
  : Plot Region of Practical Equivalence (ROPE)

## Datasets & Utilities

Example datasets and setup helpers.

- [`storks`](https://because-pkg.github.io/because/reference/storks.md)
  : Data on Stork Populations and Human Birth Rates
- [`because_family()`](https://because-pkg.github.io/because/reference/because_family.md)
  : Create a Custom Family (Distribution) for use with because
- [`because_structure()`](https://because-pkg.github.io/because/reference/because_structure.md)
  : Create a Custom Covariance Structure for use with because
- [`because_format_data()`](https://because-pkg.github.io/because/reference/because_format_data.md)
  : Format Long-Format Repeated Measures Data for Because Analysis
- [`get_family()`](https://because-pkg.github.io/because/reference/get_family.md)
  : Create a family object for a given distribution
- [`install_because_numpyro()`](https://because-pkg.github.io/because/reference/install_because_numpyro.md)
  : Install NumPyro and Python dependencies for because
