# Posterior Predictive Samples

Generic function for generating posterior predictive draws from a fitted
model. For `because` model objects, see
[`posterior_predict.because`](https://because-pkg.github.io/because/reference/posterior_predict.because.md),
which returns an `[ndraws x N_obs]` matrix of simulated response values.
These draws are the basis for
[`pp_check`](https://because-pkg.github.io/because/reference/pp_check.md).

## Usage

``` r
posterior_predict(object, ...)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Additional arguments passed to the method.

## Value

A matrix of posterior predictive draws.
