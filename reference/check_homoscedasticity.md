# Bayesian Homoscedasticity & Residual Dispersion Diagnostics

Assesses variance constancy across fitted predictions \\\hat{\mu}\\
(Gelman & Hill 2007).

Partitions fitted values into quantiles and computes the posterior
median and 95% Credible Interval of residual standard deviations in each
bin. Also computes the posterior distribution of the variance ratio
discrepancy: \$\$T\_{\text{ratio}}^{(s)} =
\frac{\text{Var}(e\_{\text{top\\bin}}^{(s)})}{\text{Var}(e\_{\text{bottom\\bin}}^{(s)})}\$\$

## Usage

``` r
# S3 method for class 'because'
check_homoscedasticity(object, resp = NULL, n_bins = 5, ndraws = 250, ...)

check_homoscedasticity(object, ...)
```

## Arguments

- object:

  A fitted model object.

- resp:

  Character string; name of the response variable. If `NULL`, takes the
  first response.

- n_bins:

  Integer; number of bins to partition fitted values into. Defaults to 5
  (quintiles).

- ndraws:

  Integer; number of posterior draws to use. Defaults to 250.

- ...:

  Additional arguments.

## Value

An object of class `"because_homoscedasticity"` containing diagnostic
plots and summary statistics.
