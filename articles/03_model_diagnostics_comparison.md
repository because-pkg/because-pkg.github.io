# Model Diagnostics and Comparison

## Introduction

After fitting a Bayesian SEM with `because`, the analysis does not end
at parameter estimation. A responsible workflow requires:

1.  **Convergence diagnostics**: Did the MCMC chains mix well?
2.  **Posterior predictive checks**: Does the model reproduce the
    observed data patterns?
3.  **Model comparison**: Which competing causal structure is best
    supported by the data?

This vignette walks through all diagnostic and comparison tools
available in `because`.

------------------------------------------------------------------------

### Setup: Fitting a Reference Model

We use a three-variable causal chain as a running example throughout.

``` r
library(because)

set.seed(42)
N <- 150

Temperature <- rnorm(N, 0, 1)    # Exogenous
NDVI        <- 0.6 * Temperature + rnorm(N, 0, 0.8)   # Mediator
Abundance   <- 0.4 * NDVI + 0.3 * Temperature + rnorm(N, 0, 1)

sim_data <- data.frame(Temperature, NDVI, Abundance)
```

``` r
# Causal model: Temperature -> NDVI -> Abundance, Temperature -> Abundance
fit <- because(
  equations = list(
    NDVI      ~ Temperature,
    Abundance ~ NDVI + Temperature
  ),
  data    = sim_data,
  dsep    = TRUE,
  WAIC    = TRUE,   # Required for model comparison
  n.iter  = 12500,
  n.burnin = 2500
)

summary(fit)
```

------------------------------------------------------------------------

### Part 1: Convergence Diagnostics

#### 1.1 The Rhat Statistic

The [`summary()`](https://rdrr.io/r/base/summary.html) output includes
**Rhat** (the Gelman-Rubin convergence diagnostic) and **n.eff** (the
effective sample size) for each parameter.

- **Rhat \< 1.05**: Excellent convergence.

- **Rhat 1.05–1.1**: Acceptable for most parameters; consider running
  more iterations.

- **Rhat \> 1.1**: Non-convergence. The chains have not mixed. Do not
  interpret these estimates.

- **n.eff \> 100**: Generally sufficient for reliable posterior
  summaries.

- **n.eff \< 100**: Estimates may be unstable. Increase `n.iter` or
  `n.thin`.

``` r
# The stored summary object has $statistics containing Rhat and n.eff
sum_stats <- fit$summary$statistics

# Identify parameters with potential convergence issues
if (!is.null(sum_stats) && "Rhat" %in% colnames(sum_stats)) {
  high_rhat <- sum_stats[sum_stats[, "Rhat"] > 1.05, , drop = FALSE]
  if (nrow(high_rhat) == 0) {
    cat("All parameters converged (Rhat < 1.05).\n")
  } else {
    print(high_rhat)
  }
} else {
  # summary() prints results to console; inspect Rhat column visually
  summary(fit)
}
```

#### 1.2 Trace Plots

Trace plots show the MCMC samples over time for each chain. A well-mixed
chain looks like a “fuzzy caterpillar”: dense, stationary, and without
long-term trends.

``` r
# Plot traces for all monitored parameters
plot(fit$samples)
```

``` r
# Focus on a specific parameter
library(coda)
plot(fit$samples[, "beta_Abundance_NDVI"])
```

**What to look for**: - Multiple chains (coloured differently)
overlapping completely. - No visible trend, drift, or “stickiness”
(chains staying at the same value for many consecutive iterations).

#### 1.3 Autocorrelation Plots

High autocorrelation within chains reduces the effective sample size.

``` r
autocorr.plot(fit$samples[, "beta_Abundance_NDVI"])
```

If autocorrelation remains high at large lags, increase `n.thin` to thin
the chains more aggressively.

------------------------------------------------------------------------

### Part 2: Posterior Predictive Checks

Posterior predictive checks (PPCs) compare the **observed data** to data
**simulated from the fitted model**. If the model is a good fit,
simulated data should resemble the real data.

`because` implements PPCs via
[`pp_check()`](https://because-pkg.github.io/because/reference/pp_check.md),
a wrapper around the `bayesplot` package.

#### 2.1 Density Overlay

``` r
# Overlay posterior predictive densities on the observed density
pp_check(fit, resp = "Abundance", type = "dens_overlay", ndraws = 50)
```

The thick line is the density of the observed data; the thin lines are
densities from 50 posterior predictive draws. They should overlap well.
Systematic mismatches (e.g., the model consistently underestimates
extreme values) indicate model misspecification.

#### 2.2 Test Statistics

``` r
# Compare the observed mean and SD to the posterior predictive distribution
pp_check(fit, resp = "Abundance", type = "stat", stat = "mean")
pp_check(fit, resp = "Abundance", type = "stat", stat = "sd")
```

The histogram shows the distribution of the test statistic (e.g., mean)
across posterior predictive datasets. The vertical line marks the
observed statistic. If the observed value falls in the tails of the
predictive distribution, the model may be misspecified for that aspect
of the data.

#### 2.3 Histogram

``` r
pp_check(fit, resp = "Abundance", type = "hist", ndraws = 10)
```

------------------------------------------------------------------------

### Part 3: Visualizing Results

#### 3.1 DAG with Path Coefficients

``` r
# Plot the causal DAG with estimated standardised path coefficients
plot_dag(fit)
```

Arrow widths are proportional to the magnitude of the standardised path
coefficients. This gives an immediate visual summary of the relative
strength of each causal effect.

#### 3.2 Coefficient Plot

``` r
# Caterpillar plot of all path coefficients with 95% credible intervals
plot_coef(fit)
```

Points are posterior means; horizontal lines span the 95% credible
interval. Coefficients whose intervals exclude zero have clear
directionality.

#### 3.3 D-Separation Results

``` r
# Visualize d-separation test results (requires dsep = TRUE at fit time)
plot_dsep(fit)
```

Each row corresponds to a conditional independence claim implied by the
causal model. Points near zero (with intervals crossing zero) indicate
that the model’s independence assumptions hold in the data.

#### 3.4 Posterior Distributions

``` r
# Compare posterior distributions for specific parameters
plot_posterior(fit, parameter = "beta_Abundance")
```

The `parameter` argument accepts partial matching (regex), so
`"beta_Abundance"` will match all path coefficients with `Abundance` as
the response variable. You can overlay multiple models:

``` r
fit_alt <- because(
  equations = list(Abundance ~ Temperature),  # Simpler model
  data = sim_data, WAIC = TRUE, quiet = TRUE
)

plot_posterior(
  model = list("Full" = fit, "Simple" = fit_alt),
  parameter = "beta_Abundance_Temperature"
)
```

------------------------------------------------------------------------

### Part 4: Information Criteria

#### 4.1 WAIC (Watanabe-Akaike Information Criterion)

WAIC is the Bayesian analogue of AIC. Lower WAIC indicates better
out-of-sample predictive performance. Enable it at fit time with
`WAIC = TRUE`:

``` r
# WAIC is reported in the summary
summary(fit)

# Or access it directly
fit$WAIC
```

The WAIC output contains:

- **WAIC**: The information criterion value (lower = better).
- **lppd**: Log pointwise predictive density (higher = better fit).
- **p_waic**: Effective number of parameters (complexity penalty).

#### 4.2 LOO-CV (Leave-One-Out Cross-Validation)

LOO-CV is generally more robust than WAIC, particularly when a small
number of observations are highly influential. It uses Pareto Smoothed
Importance Sampling (PSIS-LOO) from the `loo` package.

``` r
loo_result <- because_loo(fit)
print(loo_result)
```

##### Pareto-k Diagnostics

Each observation receives a Pareto-k value:

| k value        | Interpretation                                      |
|:---------------|:----------------------------------------------------|
| k \< 0.5       | Excellent — all estimates reliable                  |
| 0.5 ≤ k \< 0.7 | Good — estimates acceptable                         |
| 0.7 ≤ k \< 1.0 | Problematic — estimates unreliable                  |
| k ≥ 1.0        | Very problematic — refit without these observations |

``` r
# Plot Pareto-k values per observation
plot(loo_result)
```

High-k observations are influential data points that disproportionately
affect the posterior. They may indicate outliers or model
misspecification.

------------------------------------------------------------------------

### Part 5: Model Comparison

When multiple plausible causal structures exist, `because` provides
tools to rank them by predictive performance.

#### 5.1 Comparing Two Fitted Models

``` r
# Full model: Temperature -> NDVI -> Abundance, Temperature -> Abundance
fit_full <- because(
  equations = list(
    NDVI      ~ Temperature,
    Abundance ~ NDVI + Temperature
  ),
  data = sim_data, WAIC = TRUE, quiet = TRUE
)

# Simpler model: Temperature only acts through NDVI (no direct path)
fit_mediated <- because(
  equations = list(
    NDVI      ~ Temperature,
    Abundance ~ NDVI
  ),
  data = sim_data, WAIC = TRUE, quiet = TRUE
)

# Compare using WAIC
because_compare(fit_full, fit_mediated)
```

The output table ranks models by WAIC. The `deltaWAIC` column gives the
difference from the best model; `weight` is the approximate model weight
(higher = more support).

#### 5.2 Run and Compare Multiple Specifications Simultaneously

For a larger model selection exercise,
[`because_compare()`](https://because-pkg.github.io/because/reference/because_compare.md)
can run multiple model specifications in parallel and return the
comparison table directly:

``` r
model_specs <- list(
  full_mediation = list(
    equations = list(NDVI ~ Temperature, Abundance ~ NDVI + Temperature)
  ),
  partial_mediation = list(
    equations = list(NDVI ~ Temperature, Abundance ~ NDVI)
  ),
  direct_only = list(
    equations = list(Abundance ~ Temperature)
  )
)

comparison_result <- because_compare(
  model_specs = model_specs,
  data        = sim_data,
  n.cores     = 2,        # Run models in parallel
  WAIC        = TRUE,
  quiet       = TRUE
)

print(comparison_result$comparison)
```

#### 5.3 LOO-Based Comparison

For a more rigorous comparison, use LOO-CV:

``` r
loo1 <- because_loo(fit_full)
loo2 <- because_loo(fit_mediated)

# loo_compare from the loo package
library(loo)
loo_compare(loo1, loo2)
```

The `elpd_diff` column shows the difference in expected log predictive
density (positive = better). The standard error of the difference
(`se_diff`) should be considered: if `|elpd_diff| < 2 * se_diff`, the
models are approximately equivalent in predictive performance.

#### 5.4 Reusing D-Separation Tests Across Models

When comparing many alternative DAGs that share some conditional
independence claims, `because` can reuse already-computed d-separation
tests via the `reuse_models` argument, avoiding redundant JAGS runs:

``` r
fit_base <- because(
  equations = list(NDVI ~ Temperature, Abundance ~ NDVI + Temperature),
  data = sim_data, dsep = TRUE, quiet = TRUE
)

# Fit a variant that reuses any matching tests from fit_base
fit_variant <- because(
  equations = list(NDVI ~ Temperature, Abundance ~ NDVI),
  data        = sim_data,
  dsep        = TRUE,
  reuse_models = list(fit_base),  # Scan for reusable tests
  quiet       = TRUE
)
```

------------------------------------------------------------------------

### Practical Checklist

Before interpreting your results, work through this checklist:

**Rhat \< 1.05** for all parameters.

**n.eff \> 100** for all parameters of interest.

**Trace plots** show well-mixed, stationary chains.

**Posterior predictive checks** show simulated data overlapping with
observed data.

**D-separation tests** (if used) show all intervals crossing zero.

**WAIC/LOO comparison** supports the chosen causal structure over
plausible alternatives.

### References

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., &
Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). Chapman &
Hall/CRC.

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5), 1413–1432.

Shipley, B. (2016). *Cause and Correlation in Biology: A User’s Guide to
Path Analysis, Structural Equations and Causal Inference with R* (2nd
ed.). Cambridge University Press.
