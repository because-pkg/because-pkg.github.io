# Standardize Path Coefficients onto Latent Liability and Common SD Scale

Rescales parameters from a fitted `because` model or d-separation test
into standardized path coefficients (\\\beta^\*\\).

## Usage

``` r
standardize(object, ...)

# S3 method for class 'because'
standardize(
  object,
  method = c("error", "total", "mplus"),
  binary_predictors = c("raw", "sd"),
  ...
)
```

## Arguments

- object:

  A fitted `because` model object or a `because` d-separation test
  object.

- ...:

  Additional arguments passed to methods.

- method:

  Character; liability/outcome variance calculation method:

  - `"error"` (default): Disturbance error + structural/random effects
    variance (\\\sum \sigma\_{\text{random}}^2 + \pi^2/3\\ for binomial;
    \\\sum \sigma\_{\text{random}}^2 + \ln(1 + 1/\bar{Y})\\ for Poisson;
    Nakagawa et al. 2017).

  - `"total"`: Full McKelvey-Zavoina / Nakagawa total variance
    incorporating the empirical variance of the fixed linear predictor
    \\\text{Var}(\hat{\eta}\_{\text{fixed}})\\.

  - `"mplus"`: Standardizes by the standard deviation of the linear
    predictions \\s\_{xb} =
    \sqrt{\text{Var}(\hat{\eta}\_{\text{fixed}}) + \sum
    \sigma\_{\text{random}}^2}\\, matching Mplus `StdYX` for Poisson
    regression (UCLA OARC FAQ).

- binary_predictors:

  Character; scaling for binary/dummy predictors: `"raw"` (default;
  semi-standardization, keeping binary predictors on the 0-to-1 step so
  that the standardized slope represents Cohen's \\d\\ on the
  standardized outcome scale; Gelman 2008, Bollen 1989, Grace &
  Bollen 2005) or `"sd"` (full standardization, scaling by \\s_X =
  \sqrt{p(1-p)}\\; Menard 2004).

## Value

An updated object of class `c("because_standardized", "because")` with
transformed posterior samples and updated summary statistics. Fully
compatible with [`summary`](https://rdrr.io/r/base/summary.html),
[`plot_coef`](https://because-pkg.github.io/because/reference/plot_coef.md),
[`plot_rope`](https://because-pkg.github.io/because/reference/plot_rope.md),
and
[`plot_dag`](https://because-pkg.github.io/because/reference/plot_dag.md).

## Details

- **Continuous (Gaussian) Variables**: Scaled by the ratio of standard
  deviations \\\beta^\* = \beta \times (\text{SD}(X) / \text{SD}(Y))\\.
  If variables were already standardized prior to model fitting, their
  scale is preserved.

- **Binary / Categorical (Binomial/Bernoulli/Ordinal Logit)**: Raw
  coefficients sit on the log-odds scale. They are standardized *a
  posteriori* onto the continuous latent liability scale following
  McKelvey & Zavoina (1975), Menard (2004), and Nakagawa et al. (2017):
  \$\$\beta^\* = \frac{\beta \times \text{SD}(X)}{\sqrt{\sum
  \sigma\_{\text{random}}^2 + \pi^2/3}}\$\$ where \\\pi^2/3 \approx
  3.289868\\ is the theoretical disturbance variance of the standard
  logistic link.

- **Count (Poisson / Negative Binomial)**: Raw coefficients sit on the
  log-rate scale (\\\ln \lambda\\). Standardized onto the latent log /
  rate scale via:

  - `"error"` (default): \\\beta^\* = \frac{\beta \times
    \text{SD}(X)}{\sqrt{\sum \sigma\_{\text{random}}^2 + \ln(1 +
    1/\bar{Y})}}\\ following the GLMM latent distribution variance
    formulation of Nakagawa et al. (2017).

  - `"total"`: Incorporates the empirical variance of the fixed linear
    predictor \\\text{Var}(\hat{\eta}\_{\text{fixed}})\\ alongside
    random effects and disturbance variance.

  - `"mplus"`: Standardizes by the standard deviation of the linear
    predictions \\s\_{xb} =
    \sqrt{\text{Var}(\hat{\eta}\_{\text{fixed}}) + \sum
    \sigma\_{\text{random}}^2}\\, matching Mplus `StdYX` for Poisson
    regression (UCLA OARC FAQ).

- **Binary Predictors**: Slopes pointing from a binary indicator (e.g.
  `Migrant` \\\in \\0, 1\\\\) to a response are kept on their natural
  0-to-1 step when `binary_predictors = "raw"` (default), placing the
  effect on the standardized outcome scale (directly analogous to
  Cohen's \\d\\; Gelman 2008, Grace & Bollen 2005). When
  `binary_predictors = "sd"`, slopes are multiplied by \\s_X =
  \sqrt{p(1 - p)}\\ (Menard 2004).

- **d-separation Objects**: Evaluates and rescales each auxiliary
  independence claim across all d-sep sub-models, enabling a uniform
  `rope = c(-0.10, 0.10)` benchmark across all tests.

## References

McKelvey, R. D., & Zavoina, W. (1975). A statistical model for the
analysis of ordinal level dependent variables. *Journal of Mathematical
Sociology*, 4(1), 103–120.

Menard, S. (2004). Six approaches to calculating standardized logistic
regression coefficients. *The American Statistician*, 58(3), 218–223.

Gelman, A. (2008). Scaling regression inputs by dividing by two standard
deviations. *Statistics in Medicine*, 27(15), 2865–2873.

Grace, J. B., & Bollen, K. A. (2005). Interpreting the results from
multiple regression and structural equation models. *Ecology*, 86(11),
283–295.

Nakagawa, S., Johnson, P. C. D., & Schielzeth, H. (2017). The
coefficient of determination \\R^2\\ and intra-class correlation
coefficient from generalized linear mixed-effects models revisited and
expanded. *Journal of the Royal Society Interface*, 14(134), 20170213.

UCLA: Statistical Consulting Group. How does Mplus calculate the
standardized coefficients based on a Poisson model?
<https://stats.oarc.ucla.edu/mplus/faq/how-does-mplus-calculate-the-standardized-coefficients-based-on-a-poisson-model/>
