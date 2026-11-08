# Bayesian Residual Structure & Variance Partition Diagnostics

Evaluates the proportion of variance absorbed by phylogenetic or spatial
random effects in a `because` model (Hadfield 2010; de Villemereuil &
Nakagawa 2014).

Computes the posterior distribution of the **Variance Partition
Coefficient (VPC)**: \$\$\lambda\_{\text{Bayes}}^{(s)} =
\frac{\sigma\_{\text{struct}}^{2(s)}}{\sigma\_{\text{struct}}^{2(s)} +
\sigma\_{\text{res}}^{2(s)}}\$\$

## Usage

``` r
# S3 method for class 'because'
check_residual_structure(object, resp = NULL, ...)

check_residual_structure(object, ...)
```

## Arguments

- object:

  A fitted model object.

- resp:

  Character string; name of the response variable. If `NULL`, takes the
  first response.

- ...:

  Additional arguments.

## Value

An object of class `"because_residual_structure"`.
