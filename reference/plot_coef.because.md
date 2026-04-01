# plot_coef

Creates a caterpillar plot (point and whisker) of the path coefficients
from a `because` model. This provides a complementary view to the DAG,
showing the exact magnitudes and uncertainties of causal effects.

## Usage

``` r
# S3 method for class 'because'
plot_coef(
  object,
  type = "marginal",
  multinomial_probabilities = TRUE,
  color_scheme = "directional",
  ...
)
```

## Arguments

- object:

  A `because` object.

- type:

  Character; either `"marginal"` (default) or `"raw"`.

  - `"marginal"`: Shows Average Marginal Effects (AME) on the response
    scale.

  - `"raw"`: Shows the raw structural parameters (betas/rhos).

- multinomial_probabilities:

  Logical; if `TRUE` (default), expands multinomial predictors into a
  "bundle" of category-specific effects.

- color_scheme:

  Character; color scheme for significance. Options:

  - `"directional"` (default): Blue for positive significant, Red for
    negative, Grey otherwise.

  - `"sig_only"`: Black for significant, Grey for non-significant.

  - `"monochrome"`: All black.

- ...:

  Additional arguments.

## Value

A `ggplot` object.
