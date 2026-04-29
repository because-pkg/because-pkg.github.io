# Non-Gaussian Distribution Families

## Introduction

In ecology and evolutionary biology, response variables are rarely
normally distributed. Count data, binary outcomes, proportions, and
ordered categories all require appropriate distributional assumptions.
`because` natively supports a wide range of response families, all
specified through the `family` argument of
[`because()`](https://because-pkg.github.io/because/reference/because.md).

This vignette walks through each supported family with a worked
ecological example, and demonstrates how to use
[`marginal_effects()`](https://because-pkg.github.io/because/reference/marginal_effects.md)
to convert link-scale coefficients (log, logit) into interpretable
expected-value effects.

### Specifying Families

The `family` argument accepts a **named character vector**, where names
are the response variables and values are the distribution names:

``` r
library(because)

# Minimal example data
set.seed(1)
N <- 100
example_data <- data.frame(
  Y     = rpois(N, lambda = 3),
  Count = rpois(N, lambda = 5),
  Presence = rbinom(N, size = 1, prob = 0.4),
  X = rnorm(N)
)

# Single non-Gaussian response
fit <- because(
  equations = list(Y ~ X),
  data = example_data,
  family = c(Y = "poisson"),
  quiet = TRUE
)

# Multiple responses with different families
fit_mixed <- because(
  equations = list(
    Count ~ X,
    Presence ~ X
  ),
  data = example_data,
  family = c(Count = "poisson", Presence = "bernoulli"),
  quiet = TRUE
)
```

Gaussian is the default, so you only need to specify `family` for
non-Gaussian responses.

------------------------------------------------------------------------

### Example 1: Poisson — Count Data

**Use case**: Number of offspring, species counts, nest failures,
disease incidence.

**Link function**: log — `log(mu) = alpha + beta * X`

``` r
set.seed(42)
N <- 150

# Predictors
Body_Mass <- rnorm(N, mean = 50, sd = 10)  # kg
Habitat_Quality <- runif(N, 0, 1)          # 0-1 index

# True log-linear relationship
log_lambda <- 0.5 + 0.03 * Body_Mass + 1.2 * Habitat_Quality
Offspring <- rpois(N, lambda = exp(log_lambda))

count_data <- data.frame(
  Offspring = Offspring,
  Body_Mass = scale(Body_Mass)[, 1],
  Habitat_Quality = scale(Habitat_Quality)[, 1]
)
```

``` r
# Fit Poisson SEM
fit_poisson <- because(
  equations = list(Offspring ~ Body_Mass + Habitat_Quality),
  data = count_data,
  family = c(Offspring = "poisson"),
  dsep = FALSE
)

summary(fit_poisson)
```

The summary output gives `alpha_Offspring` (log-scale intercept) and
`beta_Offspring_*` (log-scale slopes). A slope of 0.5 on the log scale
means that a one-unit increase in the predictor multiplies the expected
count by $e^{0.5} \approx 1.65$.

#### Marginal Effects

To get effects on the count scale (expected number of offspring per
one-SD increase in each predictor), use
[`marginal_effects()`](https://because-pkg.github.io/because/reference/marginal_effects.md):

``` r
me_poisson <- marginal_effects(fit_poisson)
print(me_poisson)
```

The `Effect` column reports the Average Marginal Effect: the expected
change in the number of offspring per one-unit increase in the
(standardised) predictor, averaged over all observations.

------------------------------------------------------------------------

### Example 2: Bernoulli / Binomial — Presence/Absence

**Use case**: Species presence/absence, survival, reproductive success
(yes/no).

**Link function**: logit — `logit(p) = alpha + beta * X`

Use `"bernoulli"` for 0/1 binary data and `"binomial"` when you have
proportion data with known totals.

``` r
set.seed(42)
N <- 200

# Predictors
Temperature <- rnorm(N, mean = 15, sd = 5)
Precipitation <- rnorm(N, mean = 800, sd = 200)

# True logistic relationship
logit_p <- -1 + 0.15 * scale(Temperature) + 0.08 * scale(Precipitation)
Presence <- rbinom(N, size = 1, prob = 1 / (1 + exp(-logit_p)))

presence_data <- data.frame(
  Presence = Presence,
  Temperature = scale(Temperature)[, 1],
  Precipitation = scale(Precipitation)[, 1]
)
```

``` r
fit_bern <- because(
  equations = list(Presence ~ Temperature + Precipitation),
  data = presence_data,
  family = c(Presence = "bernoulli")
)

summary(fit_bern)
```

``` r
# Marginal effects: expected change in *probability* of presence per 1-SD increase
me_bern <- marginal_effects(fit_bern)
print(me_bern)
```

------------------------------------------------------------------------

### Example 3: Negative Binomial — Overdispersed Counts

**Use case**: Count data with variance that greatly exceeds the mean
(e.g., parasite loads, number of plant seeds, abundance data).

A Poisson model assumes $\text{Var}(Y) = E(Y)$. Ecological count data
often show **overdispersion** ($\text{Var}(Y) > E(Y)$), which can
inflate Type I error in Poisson models. The Negative Binomial adds a
dispersion parameter $r$ (estimated from the data) to accommodate this
extra variance.

``` r
set.seed(123)
N <- 150

Resource <- rnorm(N, 0, 1)

# Overdispersed counts: large variance relative to mean
# MASS::rnegbin: mu = exp(0.8 + 0.6*X), size = 2 (small = more overdispersed)
library(MASS)
Parasite_Load <- rnegbin(N, mu = exp(0.8 + 0.6 * Resource), theta = 2)

nb_data <- data.frame(
  Parasite_Load = Parasite_Load,
  Resource = Resource
)
```

``` r
# Compare Poisson vs. Negative Binomial
fit_pois <- because(
  equations = list(Parasite_Load ~ Resource),
  data = nb_data,
  family = c(Parasite_Load = "poisson"),
  WAIC = TRUE,
  quiet = TRUE
)

fit_nb <- because(
  equations = list(Parasite_Load ~ Resource),
  data = nb_data,
  family = c(Parasite_Load = "negbinomial"),
  WAIC = TRUE,
  quiet = TRUE
)

# Compare WAIC (lower = better fit)
because_compare(fit_pois, fit_nb)
```

The Negative Binomial model will have a substantially lower WAIC when
the data are overdispersed. The `r_Parasite_Load` parameter in the
summary captures the estimated overdispersion.

------------------------------------------------------------------------

### Example 4: Zero-Inflated Models — ZIP and ZINB

**Use case**: Count data with more zeros than expected under a standard
Poisson or Negative Binomial distribution. Common in ecological surveys
(many sites with no detections, plus genuine count variation at occupied
sites).

The Zero-Inflated Poisson (ZIP) models data as a mixture of structural
zeros and a Poisson process. The ZINB adds Negative Binomial
overdispersion.

``` r
set.seed(99)
N <- 200

Habitat <- rnorm(N, 0, 1)

# ZIP: with probability psi, the count is zero (structural zero);
# otherwise it follows a Poisson distribution.
psi <- 0.4  # 40% structural zeros
lambda <- exp(0.5 + 0.7 * Habitat)

structural_zero <- rbinom(N, 1, psi)
Abundance <- ifelse(structural_zero == 1, 0, rpois(N, lambda))

zip_data <- data.frame(
  Abundance = Abundance,
  Habitat = Habitat
)
```

``` r
fit_zip <- because(
  equations = list(Abundance ~ Habitat),
  data = zip_data,
  family = c(Abundance = "zip")
)

summary(fit_zip)
# psi_Abundance: estimated zero-inflation probability
# beta_Abundance_Habitat: effect on the Poisson component (log scale)
```

``` r
# Use ZINB when counts are also overdispersed
fit_zinb <- because(
  equations = list(Abundance ~ Habitat),
  data = zip_data,
  family = c(Abundance = "zinb")
)

summary(fit_zinb)
```

------------------------------------------------------------------------

### Example 5: Ordinal — Ordered Categorical Responses

**Use case**: Behaviour scores, habitat quality rankings, dominance
rank, pain scales, Likert-scale data.

Ordinal models estimate a set of **cutpoints** (thresholds) that divide
a latent continuous scale into categories. The proportional odds
assumption means a single slope parameter describes the effect on all
adjacent-category boundaries.

``` r
set.seed(7)
N <- 200

# Predictor: years of experience
Dominance_Rank <- rnorm(N, 0, 1)

# Underlying latent variable
latent <- 0 + 1.5 * Dominance_Rank + rnorm(N, 0, 1)

# 4 behaviour categories: Subordinate (1), Low (2), Medium (3), Dominant (4)
Behaviour <- cut(latent,
  breaks = c(-Inf, -1, 0, 1, Inf),
  labels = FALSE
)

ord_data <- data.frame(
  Behaviour = Behaviour,
  Dominance_Rank = Dominance_Rank
)
```

``` r
fit_ord <- because(
  equations = list(Behaviour ~ Dominance_Rank),
  data = ord_data,
  family = c(Behaviour = "ordinal")
)

summary(fit_ord)
# cutpoint_Behaviour[k]: thresholds on the latent scale separating categories k and k+1
# beta_Behaviour_Dominance_Rank: effect on the latent scale (positive = higher categories)
```

#### The `expand_ordered` argument

By default, `because` uses a proportional-odds parameterisation with a
single slope for all category boundaries. If you want
**category-specific slopes** (relaxing the proportional odds
assumption), set `expand_ordered = TRUE`:

``` r
fit_ord_expanded <- because(
  equations = list(Behaviour ~ Dominance_Rank),
  data = ord_data,
  family = c(Behaviour = "ordinal"),
  expand_ordered = TRUE  # Separate slope per adjacent-category boundary
)

summary(fit_ord_expanded)
```

This produces separate `beta_Behaviour_Dominance_Rank[k]` parameters for
each cutpoint boundary.

#### Marginal Effects for Ordinal Responses

``` r
# Expected category (1–4) per 1-unit increase in Dominance_Rank
me_ord <- marginal_effects(fit_ord)
print(me_ord)
```

------------------------------------------------------------------------

### Summary: Choosing a Family

| Data Type                     | Family                 | Link               | Key Parameter             |
|:------------------------------|:-----------------------|:-------------------|:--------------------------|
| Continuous, symmetric         | `"gaussian"` (default) | identity           | `sigma_e_*`               |
| Counts, mean ≈ variance       | `"poisson"`            | log                | —                         |
| Counts, variance \>\> mean    | `"negbinomial"`        | log                | `r_*`                     |
| Counts with excess zeros      | `"zip"`                | log (count part)   | `psi_*`                   |
| Counts, overdispersed + zeros | `"zinb"`               | log (count part)   | `r_*`, `psi_*`            |
| Binary (0/1)                  | `"bernoulli"`          | logit              | —                         |
| Proportions with totals       | `"binomial"`           | logit              | —                         |
| Ordered categories            | `"ordinal"`            | logit (cumulative) | `cutpoint_*[k]`           |
| Unordered categories          | `"multinomial"`        | log (softmax)      | `alpha_*[k]`, `beta_*[k]` |

### A Note on Convergence for Non-Gaussian Models

Non-Gaussian link functions (log, logit) create a curved relationship
between the linear predictor and the likelihood, which can make MCMC
mixing slower. Recommendations:

1.  **Standardize predictors** before fitting:
    `data = as.data.frame(scale(mydata))`.
2.  **Use more iterations**: `n.iter = 25000` or more.
3.  **Check Rhat**: all parameters should have Rhat \< 1.05.
4.  **Use the NIMBLE engine** for better mixing on non-Gaussian models:
    `engine = "nimble"` (see the [NIMBLE
    vignette](https://because-pkg.github.io/because/articles/10_nimble_engine.md)).

### References

Agresti, A. (2013). *Categorical Data Analysis* (3rd ed.). Wiley.

Martin, T. G., Wintle, B. A., Rhodes, J. R., Kuhnert, P. M., Field, S.
A., Low-Choy, S. J., … Possingham, H. P. (2005). Zero tolerance ecology:
improving ecological inference by modelling the source of zero
observations. *Ecology Letters*, 8(11), 1235–1246.
<https://doi.org/10.1111/j.1461-0248.2005.00826.x>

McCullagh, P., & Nelder, J. A. (1989). *Generalized Linear Models* (2nd
ed.). Chapman & Hall.

Ver Hoef, J. M., & Boveng, P. L. (2007). Quasi-Poisson vs. negative
binomial regression: How should we model overdispersed count data?
*Ecology*, 88(11), 2766–2772. <https://doi.org/10.1890/07-0043.1>

Zuur, A. F., Ieno, E. N., Walker, N. J., Saveliev, A. A., & Smith, G. M.
(2009). *Mixed Effects Models and Extensions in Ecology with R*.
Springer.
