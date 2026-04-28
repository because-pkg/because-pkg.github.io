# Posterior Predictive Checks

Generic function for posterior predictive checks of a fitted model. For
`because` model objects, see
[`pp_check.because`](https://because-pkg.github.io/because/reference/pp_check.because.md),
which wraps `bayesplot` functions (density overlay, histogram, test
statistics) and supports conditional or marginal prediction via
`re_formula`.

## Usage

``` r
pp_check(object, ...)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Additional arguments passed to the method.

## Value

A ggplot object.
