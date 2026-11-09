# Bayesian Normality Check & Q-Q Envelope

Evaluates distributional and normality assumptions using **Randomized
Quantile Residuals** (Dunn & Smyth 1996) and **Posterior Predictive
Credible Envelopes** (Landwehr et al. 1984).

Generates an empirical simulation envelope (50%, 80%, and 95% bands)
directly from posterior predictive replications \\y^{\text{rep}}\\, and
evaluates sample skewness and kurtosis test quantities against their
posterior predictive distributions.

## Usage

``` r
# S3 method for class 'because'
check_normality(object, resp = NULL, ndraws = 250, ...)

check_normality(object, ...)
```

## Arguments

- object:

  A fitted model object.

- resp:

  Character vector; name(s) of response variable(s) to evaluate. If
  `NULL` (default), evaluates all endogenous response variables in the
  model.

- ndraws:

  Integer; number of posterior predictive draws to simulate for the
  envelope. Defaults to 250.

- ...:

  Additional arguments.

## Value

An object of class `"because_normality"` containing the ggplot object,
envelope data, and skewness/kurtosis posterior tail probabilities.
