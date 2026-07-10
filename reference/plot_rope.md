# Plot Region of Practical Equivalence (ROPE)

Generic for plotting posterior density distributions and visualizing
their overlap with a Region of Practical Equivalence (ROPE).

## Usage

``` r
plot_rope(object, ...)

# S3 method for class 'because'
plot_rope(object, rope = c(-0.1, 0.1), parameters = NULL, ...)
```

## Arguments

- object:

  A `because` object.

- ...:

  Additional arguments passed to methods.

- rope:

  A numeric vector of length 2 specifying the ROPE limits. Default is
  `c(-0.1, 0.1)`.

- parameters:

  Optional character vector to filter which parameters to plot.

## Value

A `ggplot` object.
