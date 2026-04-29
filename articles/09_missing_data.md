# Bayesian Missing Data Imputation

## Introduction

Missing data are ubiquitous in ecology and evolutionary biology. Traits
may not be measurable for all individuals, observations may be lost, or
certain variables may be systematically unavailable for a subset of
observations. How missing data are handled can profoundly affect the
quality of causal inference.

### Missing Data Mechanisms

Before choosing an imputation strategy, it is important to understand
*why* data are missing (Rubin, 1976; Little & Rubin, 2002):

- **MCAR (Missing Completely At Random)**: Missingness is unrelated to
  any variable in the study, observed or unobserved. The missing rows
  are a random subsample. List-wise deletion is unbiased under MCAR but
  wasteful.
- **MAR (Missing At Random)**: Missingness depends on observed variables
  but not on the unobserved values themselves. For example, heavier
  animals may be harder to weigh, but once body length is accounted for,
  the probability of a missing mass value does not depend on the true
  mass. Model-based imputation (including Bayesian imputation) is valid
  under MAR.
- **MNAR (Missing Not At Random)**: Missingness depends on the
  unobserved values themselves (e.g., extremely light animals die before
  weighing). No standard imputation method can fully correct for MNAR
  without external information about the missing-data mechanism.

`because` handles missing data through **full Bayesian imputation**:
missing values are treated as unknown parameters and estimated
simultaneously with all other model parameters during MCMC sampling.
This approach has important advantages over common alternatives such as
list-wise deletion or two-step multiple imputation (van Buuren &
Groothuis-Oudshoorn, 2011):

1.  **No loss of power**: All observations contribute to parameter
    estimation, even rows with missing predictors.
2.  **Logical consistency**: Because imputation happens inside the
    model, imputed values always respect the causal structure (e.g.,
    deterministic relationships between variables are preserved in every
    MCMC draw).
3.  **Valid under MAR**: Like other model-based approaches, Bayesian
    joint-model imputation provides unbiased estimates when data are
    MCAR or MAR. Under MNAR, bias can remain; consider adding covariates
    that predict missingness to reduce (but not eliminate) this bias
    (Nakagawa & Freckleton, 2008).

This vignette demonstrates how to use `because`’s built-in imputation
and the
[`extract_imputed()`](https://because-pkg.github.io/because/reference/extract_imputed.md)
function to recover the posterior distribution of imputed values.

------------------------------------------------------------------------

### Why Bayesian Imputation Within the SEM?

Consider a dataset where `IsMature` is defined as `Age > 2`. A standard
two-step approach (impute first, then model) would:

1.  Impute `Age` independently using regression on other variables.
2.  Impute `IsMature` independently based on hormone levels or body
    size.
3.  **Result**: A marmot aged 6 that is classified as Immature —
    physically impossible.

With `because`, `IsMature` is a **deterministic node** defined by `Age`.
The model always calculates `IsMature` as a function of the imputed
`Age` in every MCMC step, guaranteeing logical coherence throughout the
posterior. See the [Deterministic Nodes
vignette](https://because-pkg.github.io/because/articles/06_deterministic_nodes.md)
for details.

------------------------------------------------------------------------

### Example 1: Missing Predictors

We simulate an ecological dataset where individual body mass (`Mass`)
affects annual reproductive output (`Offspring`), but mass is missing
for 20% of individuals.

#### Simulate Data with Missing Values

``` r
library(because)

set.seed(42)
N <- 100

# Predictors (fully observed)
Age <- runif(N, 1, 10)
Site_Quality <- rnorm(N, 0, 1)

# Body mass (partially observed)
Mass <- 5 + 0.8 * Age + 1.2 * Site_Quality + rnorm(N, 0, 1.5)

# Outcome: annual reproductive output
Offspring <- rpois(N, lambda = exp(0.3 + 0.15 * Mass))

# Introduce ~20% missingness in Mass
miss_idx <- sample(1:N, size = round(0.2 * N))
Mass_obs <- Mass
Mass_obs[miss_idx] <- NA

eco_data <- data.frame(
  Offspring = Offspring,
  Mass = Mass_obs,        # Has NAs
  Age = Age,
  Site_Quality = Site_Quality
)

# How many are missing?
cat("Missing Mass values:", sum(is.na(eco_data$Mass)), "\n")
```

#### Fit the Model

Missing values in the data frame are automatically detected and imputed
during model fitting. **No special syntax is required** — just include
the column with NAs as usual:

``` r
equations <- list(
  Mass     ~ Age + Site_Quality,  # Mass depends on Age and Site quality
  Offspring ~ Mass                # Offspring depends on Mass (count outcome)
)

# Run with monitor = "all" to record imputed values in the posterior
fit_imp <- because(
  equations = equations,
  data      = eco_data,
  family    = c(Offspring = "poisson"),
  monitor   = "all",     # Essential: saves imputed node posteriors
  n.iter    = 15000,
  n.burnin  = 3000
)

summary(fit_imp)
```

> **Note:** `monitor = "all"` is required to record the posterior
> distributions of imputed values. With the default
> `monitor = "interpretable"`, only regression coefficients and variance
> components are saved, and
> [`extract_imputed()`](https://because-pkg.github.io/because/reference/extract_imputed.md)
> will find no imputed nodes.

#### Extract Imputed Values

``` r
# Extract posterior summaries for all imputed Mass values
imputed <- extract_imputed(fit_imp)

# View the first few rows
head(imputed)
#   Variable RowIndex     Mean       SD      Q2.5      Q50     Q97.5
# 1     Mass        3   7.234    0.843     5.595    7.218     8.901
# 2     Mass       11   6.119    0.931     4.336    6.098     7.998
# ...
```

The output data frame contains:

- **Variable**: The variable with the missing value (`Mass`).
- **RowIndex**: The row number in the original data.
- **Mean, SD**: Posterior mean and standard deviation of the imputed
  value.
- **Q2.5, Q50, Q97.5**: 2.5%, 50%, and 97.5% credible interval bounds.

#### Visualize Imputed Values

We can plot the imputed values with their credible intervals alongside
the observed data to assess imputation quality:

``` r
# Add a column indicating which rows were imputed
eco_data$imputed <- is.na(eco_data$Mass)

# Merge imputed summaries back onto the data
# (imputed$RowIndex gives row positions in the original data)
plot_data <- eco_data
plot_data$Mass_mean <- plot_data$Mass
for (i in seq_len(nrow(imputed))) {
  idx <- imputed$RowIndex[i]
  plot_data$Mass_mean[idx] <- imputed$Mean[i]
}

library(ggplot2)
ggplot(plot_data, aes(x = Age, y = Mass_mean, colour = imputed)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "tomato"),
    labels = c("Observed", "Imputed")
  ) +
  labs(
    x = "Age", y = "Body Mass",
    title = "Observed and Bayesian-Imputed Body Mass",
    colour = ""
  ) +
  theme_minimal()
```

Imputed values (red) should fall within the plausible range given the
model — not at the boundary of observed values or at implausible
extremes.

------------------------------------------------------------------------

### Example 2: Missing Response Variable

`because` can also impute missing values on the **response** side. This
is useful when you want to infer the likely values of an unmeasured
outcome for particular individuals.

``` r
set.seed(7)

# Some Offspring values are missing (field data collection gaps)
eco_data2 <- eco_data
eco_data2$Offspring[sample(1:N, 10)] <- NA

fit_imp2 <- because(
  equations = list(
    Mass     ~ Age + Site_Quality,
    Offspring ~ Mass
  ),
  data    = eco_data2,
  family  = c(Offspring = "poisson"),
  monitor = "all",
  n.iter  = 15000,
  n.burnin = 3000,
  quiet   = TRUE
)

imputed2 <- extract_imputed(fit_imp2)
print(imputed2)
```

------------------------------------------------------------------------

### Example 3: Missing Data and Deterministic Nodes

As mentioned in the introduction, Bayesian imputation preserves
deterministic constraints. Here we demonstrate this with the marmot
example from the [Deterministic Nodes
vignette](https://because-pkg.github.io/because/articles/06_deterministic_nodes.md):

``` r
set.seed(99)
N <- 80

Age <- runif(N, 0, 8)
IsMature <- as.numeric(Age >= 2)  # Deterministic: 1 if age >= 2
Mating_Success <- 2.5 * IsMature + rnorm(N, 0, 0.8)

# Introduce missingness in Age
Age_obs <- Age
Age_obs[sample(1:N, 15)] <- NA

marmot_data <- data.frame(
  Mating_Success = Mating_Success,
  Age = Age_obs
)

# The formula I(Age >= 2) defines IsMature as a deterministic node.
# For every MCMC draw, if Age is imputed as 5, IsMature is automatically set to 1.
# There is no inconsistency between Age and IsMature in any posterior sample.
fit_marmot <- because(
  equations = list(Mating_Success ~ I(Age >= 2)),
  data      = marmot_data,
  monitor   = "all",
  n.iter    = 15000,
  n.burnin  = 3000,
  quiet     = TRUE
)

summary(fit_marmot)
```

------------------------------------------------------------------------

### Comparison: Imputation vs. List-wise Deletion

To illustrate the cost of discarding incomplete observations, we fit the
same model using only complete cases:

``` r
# List-wise deletion
complete_data <- eco_data[!is.na(eco_data$Mass), ]
cat("Remaining observations:", nrow(complete_data), "of", N, "\n")

fit_listwise <- because(
  equations = list(
    Mass     ~ Age + Site_Quality,
    Offspring ~ Mass
  ),
  data   = complete_data,
  family = c(Offspring = "poisson"),
  quiet  = TRUE
)

# Compare posterior widths
summary(fit_imp)      # Uses all N observations (some imputed)
summary(fit_listwise) # Uses only complete cases
```

The imputation model typically produces **narrower credible intervals**
for the path coefficients because it uses more information. Under MAR
(the standard assumption), it also avoids the bias from non-random
missingness that would arise with list-wise deletion on predictors.
Under MNAR, however, both list-wise deletion and Bayesian imputation can
be biased; the Bayesian approach is still preferred because it retains
power and allows sensitivity analyses.

------------------------------------------------------------------------

### Practical Guidelines

| Situation                                        | Recommendation                                                                                                                                    |
|:-------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------|
| Predictors have NAs                              | Include them as-is; `because` auto-imputes. Use `monitor = "all"` to extract posteriors.                                                          |
| Response variable has NAs                        | Same — NAs in any column are handled automatically.                                                                                               |
| Missingness is MCAR or MAR                       | Bayesian joint imputation provides valid, efficient inference.                                                                                    |
| Missingness is informative (MNAR)                | Consider adding covariates that predict missingness to the model; conduct sensitivity analyses. No standard method fully removes MNAR bias.       |
| Deterministic relationships among variables      | Use [Deterministic Nodes](https://because-pkg.github.io/because/articles/06_deterministic_nodes.md); constraints are preserved during imputation. |
| Very high missingness (\>50%) in a key predictor | Inspect model fit carefully; imputation becomes unreliable without sufficient observed anchor points.                                             |

### References

Little, R. J. A., & Rubin, D. B. (2002). *Statistical Analysis with
Missing Data* (2nd ed.). Wiley.

Nakagawa, S., & Freckleton, R. P. (2008). Missing inaction: the dangers
of ignoring missing data. *Trends in Ecology & Evolution*, 23(11),
592–596. <https://doi.org/10.1016/j.tree.2008.06.014>

Rubin, D. B. (1976). Inference and missing data. *Biometrika*, 63(3),
581–592. <https://doi.org/10.1093/biomet/63.3.581>

van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate
imputation by chained equations in R. *Journal of Statistical Software*,
45(3), 1–67. <https://doi.org/10.18637/jss.v045.i03>
