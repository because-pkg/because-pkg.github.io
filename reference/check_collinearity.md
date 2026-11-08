# Bayesian Collinearity & Variance Inflation Diagnostics

Evaluates collinearity among predictors in each equation of a fitted
`because` model.

Computes the **Posterior Parameter Correlation Matrix**
\\\text{Cor}(\boldsymbol{\beta} \mid y)\\ and **Leamer's Bayesian
Variance Inflation Factor (BVIF)**: \$\$\text{BVIF}\_j =
\[\boldsymbol{\Sigma}\_{\boldsymbol{\beta}}\]\_{jj} \cdot
\[\boldsymbol{\Sigma}\_{\boldsymbol{\beta}}^{-1}\]\_{jj} = \frac{1}{1 -
R^2\_{\beta_j \mid \boldsymbol{\beta}\_{-j}}}\$\$ where
\\\boldsymbol{\Sigma}\_{\boldsymbol{\beta}}\\ is the posterior
covariance matrix of the regression parameters.

The square root of BVIF gives the **SE Multiplier**, quantifying the
expansion in posterior uncertainty attributable to shared informational
content between predictors.

## Usage

``` r
# S3 method for class 'because'
check_collinearity(object, equation = NULL, ...)

bvif(object, ...)

check_collinearity(object, ...)
```

## Arguments

- object:

  A fitted model object.

- equation:

  Optional character string or integer specifying a single equation to
  evaluate. If `NULL` (default), evaluates all equations with 2 or more
  predictors.

- ...:

  Additional arguments.

## Value

An object of class `"because_collinearity"`, containing summary tables
and correlation matrices.

## Details

**Causal Notice:** In a causal DAG, elevated variance inflation between
confounders or mediators is an expected consequence of network topology.
Variables required by the DAG's adjustment set should **not** be dropped
based on BVIF, as doing so introduces confounding bias.
