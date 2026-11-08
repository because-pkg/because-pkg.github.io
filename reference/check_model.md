# Unified Diagnostic Dashboard for because Models

Generates a comprehensive multi-panel visual diagnostic summary for a
`because` model equation.

Panels included:

1.  **Posterior Predictive Check** (`pp_check`): Observed density vs.
    replicated densities.

2.  **Bayesian Q-Q Plot** (`check_normality`): Randomized quantile
    residuals with credible simulation envelope.

3.  **Residuals vs. Fitted** (`check_homoscedasticity`): Pointwise
    residual spread.

4.  **Binned Residual Variance** (`check_homoscedasticity`): Variance
    dispersion across predictions.

5.  **Collinearity & Variance Inflation** (`check_collinearity`):
    Predictor BVIF and uncertainty inflation.

## Usage

``` r
# S3 method for class 'because'
check_model(object, resp = NULL, ndraws = 250, ...)

check_model(object, ...)
```

## Arguments

- object:

  A fitted model object.

- resp:

  Character string; name of the response variable. If `NULL`, takes the
  first response.

- ndraws:

  Integer; number of draws to use for predictive checks. Defaults to
  250.

- ...:

  Additional arguments.

## Value

A list containing all diagnostic sub-objects and plots, printed
invisibly.
