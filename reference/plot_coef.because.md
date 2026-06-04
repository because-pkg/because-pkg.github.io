# Plot Path Coefficients (Caterpillar Plot)

Creates a 'caterpillar plot' (point and whisker) of all path
coefficients from a `because` model. This visualization provides a
precise statistical complement to the
[`plot_dag()`](https://because-pkg.github.io/because/reference/plot_dag.md)
overview, allowing for a side-by-side comparison of effect sizes and
their Bayesian credibility intervals.

## Usage

``` r
# S3 method for class 'because'
plot_coef(
  object,
  type = "raw",
  multinomial_probabilities = TRUE,
  color_scheme = "sig_only",
  ...
)
```

## Arguments

- object:

  A `because` object.

- type:

  Character; either `"raw"` (default) or `"marginal"`.

  - `"marginal"`: Shows **Average Marginal Effects (AME)**. For
    categorical predictors, this represents the average shift in the
    outcome (e.g. probability or counts) associated with a one-category
    change. This is the recommended scale for comparing cross-model
    impacts.

  - `"raw"` (default): Shows the raw structural parameters (betas/rhos)
    from the JAGS model. Useful for model diagnostics but harder to
    interpret on the original data scale.

- multinomial_probabilities:

  Logical; if `TRUE` (default), expands multinomial predictors into a
  "bundle" of category-specific effects, matching the arcs in the DAG.

- color_scheme:

  Character; color scheme for significance. Options:

  - `"sig_only"` (default): Discrete Black/Grey scheme. Significant
    paths (where 95% CI excludes zero) are Black; non-significant are
    Light Grey.

  - `"directional"`: Switched to a directional Red/Blue/Grey scheme.

  - `"monochrome"`: All effects are rendered in Black regardless of
    significance.

- ...:

  Additional arguments.

## Value

A `ggplot` object. Use standard `ggplot2` functions like `+ ggtitle()`
or `+ theme()` to further customize the output.

## Details

The function automatically sorts coefficients by the **Response
Variable**, effectively grouping all predictors for a given outcome
together. This hierarchy makes it intuitive to read 'down' the plot to
see what factors contribute most to a specific part of your causal
system.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit a model
fit <- because(list(Y ~ X + Z, X ~ Z), family = c(Y="binomial", X="gaussian"), data = dat)

# Plot marginal effects for intuitive interpretation
plot_coef(fit, type = "marginal")

# Customizing the plot
library(ggplot2)
plot_coef(fit) + labs(title = "Causal Influence on Wildlife Tolerance")
} # }
```
