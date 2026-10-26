# Plot Region of Practical Equivalence (ROPE)

Generic for plotting posterior density distributions and visualizing
their overlap with a Region of Practical Equivalence (ROPE).

## Usage

``` r
plot_rope(object, ...)

# S3 method for class 'because'
plot_rope(object, rope = "default", parameters = NULL, ...)
```

## Arguments

- object:

  A `because` object.

- ...:

  Additional arguments passed to methods.

- rope:

  Optional ROPE specification. Can be `"default"` (uses principled
  family defaults: `c(-0.10, 0.10)` for Gaussian, `c(-0.05, 0.05)` for
  Poisson/NegBin, `c(-0.18, 0.18)` for Binomial/Ordinal/Multinomial), a
  numeric vector of length 2 (e.g., `c(-0.1, 0.1)` applied globally), or
  a named list specifying custom intervals by family (e.g.,
  `list(poisson = c(-0.05, 0.05))`), by response variable, or by
  parameter name.

- parameters:

  Optional character vector to filter which parameters to plot, or
  `"all"` to plot all monitored parameters. If `NULL` (default), only
  structural path coefficients (`beta_*`) and correlations (`rho_*`) are
  plotted for standard models, excluding nuisance parameters (intercepts
  and variances).

## Value

A `ggplot` object.

## Details

The Region of Practical Equivalence (ROPE) formalizes Bayesian
equivalence testing (Kruschke 2018; Makowski et al. 2019). When
evaluating causal paths or conditional independence claims in
d-separation, failure to reject a sharp null (\\\beta = 0\\) does not
provide evidence of absence. ROPE defines an explicit negligible
boundary \\\[-\delta, +\delta\]\\ around zero:

- **Gaussian**: The default interval is `c(-0.10, 0.10)`, representing
  half of Cohen's negligible effect size benchmark (\\\pm 0.1 \times
  \text{SD}\\) for standardized continuous traits (Kruschke 2018).

- **Poisson & Negative Binomial**: Coefficients sit on the log-rate
  scale (\\\log \text{IRR}\\). The default interval `c(-0.05, 0.05)`
  represents an approximate \\\pm 5\\\\ relative change in expected
  count (\\\text{IRR} \in \[\exp(-0.05), \exp(0.05)\] \approx \[0.95,
  1.05\]\\).

- **Binomial, Ordinal, & Multinomial**: In logistic models, the binary
  outcome is modeled via an underlying continuous latent variable with
  standard logistic error variance \\\sigma^2 = \pi^2 / 3\\ and standard
  deviation \\\sigma = \pi / \sqrt{3} \approx 1.814\\. Applying
  Kruschke's \\0.1 \times \sigma\\ rule to this latent scale yields
  \\\pm 0.1 \times (\pi / \sqrt{3}) \approx \pm 0.181\\ (Gelman et al.
  2008; Makowski et al. 2019). This corresponds to an Odds Ratio window
  of \\\text{OR} \in \[\exp(-0.181), \exp(0.181)\] \approx \[0.83,
  1.20\]\\.

## References

Kruschke, J. K. (2018). Rejecting or accepting parameter values in
Bayesian estimation. *Advances in Methods and Practices in Psychological
Science*, 1(2), 270–280.
[doi:10.1177/2515245918771304](https://doi.org/10.1177/2515245918771304)

Makowski, D., Ben-Shachar, M. S., & Lüdecke, D. (2019). bayestestR:
Describing effects and their uncertainty, existence and significance
within the Bayesian framework. *Journal of Open Source Software*, 4(40),
1541. [doi:10.21105/joss.01541](https://doi.org/10.21105/joss.01541)

Gelman, A., Jakulin, A., Pittau, M. G., & Su, Y. S. (2008). A weakly
informative default prior distribution for logistic and other regression
models. *The Annals of Applied Statistics*, 2(4), 1360–1383.
[doi:10.1214/08-AOAS191](https://doi.org/10.1214/08-AOAS191)
