# Plot ROPE for Because Model

Plot ROPE for Because Model

## Usage

``` r
# S3 method for class 'because'
plot_rope(object, rope = c(-0.1, 0.1), parameters = NULL, ...)
```

## Arguments

- object:

  A `because` object.

- rope:

  A numeric vector of length 2 specifying the ROPE limits. Default is
  `c(-0.1, 0.1)`.

- parameters:

  Optional character vector to filter which parameters to plot.

- ...:

  Additional arguments passed to methods.

## Value

A `ggplot` object.
