# Summary for Because Model

Summarizes the output of a Because model run.

## Usage

``` r
# S3 method for class 'because'
summary(
  object,
  show_internal = FALSE,
  show_nodes = FALSE,
  show_random = FALSE,
  prob = 0.95,
  rope = NULL,
  ...
)
```

## Arguments

- object:

  A `because` object.

- show_internal:

  Logical. If `TRUE`, shows internal parameters created for
  deterministically defined nodes (e.g., `beta_..._det_...`). Defaults
  to `FALSE`.

- show_nodes:

  Logical. If `TRUE`, shows latent node values (e.g., `Age[1]`).
  Defaults to `FALSE` to prevent clutter when `monitor="all"`.

- show_random:

  Logical. If `TRUE`, shows random effect estimates (e.g., `u_...`).
  Defaults to `FALSE`.

- prob:

  Numeric. The credibility interval probability mass. Defaults to 0.95.

- rope:

  Optional ROPE specification. Can be `NULL` (default), `TRUE` (uses
  principled family defaults: `c(-0.10, 0.10)` for Gaussian,
  `c(-0.05, 0.05)` for Poisson/NegBin, `c(-0.18, 0.18)` for
  Binomial/Ordinal/Multinomial), a numeric vector of length 2 (e.g.,
  `c(-0.1, 0.1)` applied globally), or a named list specifying custom
  intervals by family (e.g., `list(poisson = c(-0.05, 0.05))`), by
  response variable, or by parameter name. If provided, calculates the
  proportion of posterior samples falling inside the ROPE (`P_ROPE`).

- ...:

  Additional arguments passed to
  [`summary.mcmc`](https://rdrr.io/pkg/coda/man/summary.mcmc.html).

## Value

A summary object containing statistics for the monitored parameters. If
`dsep = TRUE` was used in `because`, the summary focuses on the
conditional independence tests.

## Details

When `rope` is specified (or `rope = TRUE`), the function calculates
`P_ROPE`, the proportion of posterior draws falling within the Region of
Practical Equivalence. By default, principled family-specific defaults
are applied:

- **Gaussian**: `c(-0.10, 0.10)`, representing \\\pm 0.10 \times
  \text{SD}\\ for standardized continuous traits (Kruschke 2018).

- **Poisson & Negative Binomial**: `c(-0.05, 0.05)`, representing an
  approximate \\\pm 5\\\\ relative change in expected count
  (\\\text{IRR} \in \[0.95, 1.05\]\\).

- **Binomial, Ordinal, & Multinomial**: `c(-0.18, 0.18)`, representing
  \\\pm 0.10 \times (\pi / \sqrt{3})\\ on the latent logistic scale
  (Odds Ratio \\\text{OR} \in \[0.83, 1.20\]\\; Gelman et al. 2008;
  Makowski et al. 2019).

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

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- because(list(Y ~ X), data = my_data)

# Standard summary
summary(fit)

# Show latent states and random effects
summary(fit, show_nodes = TRUE, show_random = TRUE)
} # }
```
