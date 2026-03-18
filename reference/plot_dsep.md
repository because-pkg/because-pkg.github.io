# plot_dsep

Creates a caterpillar plot (point and whisker) of the regression
coefficients from all d-separation tests. A horizontal red line at zero
helps visually assess which independence claims are fulfilled (95% CI
includes zero) or violated (95% CI excludes zero).

## Usage

``` r
# S3 method for class 'because'
plot_dsep(object, ...)

plot_dsep(object, ...)
```

## Arguments

- object:

  A `because` object fitted with `dsep = TRUE`.

- ...:

  Additional arguments.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot results for a fitted model
plot_dsep(fit)
} # }
```
