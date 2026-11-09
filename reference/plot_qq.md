# Plot Bayesian Residual Q-Q Plots with Credible Simulation Envelopes

Evaluates distributional and normality assumptions for endogenous
variables in a `because` model using **Randomized Quantile Residuals**
(Dunn & Smyth 1996) and **Posterior Predictive Credible Envelopes**
(Landwehr et al. 1984).

## Usage

``` r
# S3 method for class 'because'
plot_qq(object, resp = NULL, ndraws = 250, ...)

plot_qq(object, ...)
```

## Arguments

- object:

  A fitted `because` model object.

- resp:

  Character vector of response variable names to evaluate. If `NULL`
  (the default), all endogenous response variables in the model are
  plotted and faceted into a multi-panel display. Provide a single
  variable name (e.g. `resp = "Brain"`) to plot a single variable, or a
  character vector to select a subset of responses.

- ndraws:

  Integer; number of posterior predictive draws to simulate for the
  envelope. Defaults to 250.

- ...:

  Additional arguments passed to methods.

## Value

A `ggplot2` plot object.

## Details

Standard Q-Q plots for non-Gaussian or structured models can be
misleading because discrete or constrained responses do not follow
standard normal error distributions. `plot_qq` resolves this by
computing **Randomized Quantile Residuals** via the Probability Integral
Transform (PIT) against simulated model replications
(\\y^{\text{rep}}\\), which follow an exact standard normal distribution
under the correct model (Dunn & Smyth 1996; Hartig 2022).

In addition, an empirical simulation envelope is constructed by sorting
replicated residuals across posterior predictive draws (Landwehr et al.
1984):

- The **darker blue central band** represents the 50% posterior
  predictive envelope (interquartile range).

- The **lighter blue outer band** represents the 95% posterior
  predictive envelope.

Residual points falling systematically outside the 95% envelope indicate
model misspecification, unmodeled skewness, heavy tails, or
zero-inflation.

When the model contains multiple equations (as in structural equation
models), calling `plot_qq(fit)` automatically facets across all
endogenous response variables.

## References

Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics*, 5(3), 236-244.

Landwehr, J. M., Pregibon, D., & Shoemaker, A. C. (1984). Graphical
methods for assessing logistic regression models. *Journal of the
American Statistical Association*, 79(385), 61-71.

Hartig, F. (2022). DHARMa: Residual Diagnostics for Hierarchical
(Multi-Level / Mixed) Regression Models.

Gelman, A., Meng, X. L., & Stern, H. (1996). Posterior predictive
assessment of model fitness via realized discrepancies. *Statistica
Sinica*, 733-760.

## See also

[`check_normality`](https://because-pkg.github.io/because/reference/check_normality.md),
[`check_model`](https://because-pkg.github.io/because/reference/check_model.md),
[`residuals.because`](https://because-pkg.github.io/because/reference/residuals.because.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot Q-Q envelopes for all equations in the SEM
plot_qq(fit_model)

# Plot a specific response variable
plot_qq(fit_model, resp = "Brain")
} # }
```
