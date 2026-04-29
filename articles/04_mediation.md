# Mediation Analysis

## Introduction

Causal inference often requires investigating *how* an effect occurs.
Does $X$ affect $Y$ directly, or does it work through a mediator $M$?

The `because` package provides a fully automated Bayesian mediation
analysis tool,
[`because_mediation()`](https://because-pkg.github.io/because/reference/because_mediation.md),
which decomposes the Total Effect of an exposure on an outcome into: 1.
**Direct Effect**: The effect of $\left. X\rightarrow Y \right.$ not
mediated by other variables in the graph. 2. **Indirect Effect(s)**: The
effect propagated through intermediate variables
($\left. X\rightarrow M\rightarrow Y \right.$).

This is calculated by multiplying the posterior distributions of
coefficients along each path, preserving full uncertainty
quantification.

## Example: Ecological Mediation (Elevation Gradient)

In this example, we investigate how **Elevation** affects **Plant
Abundance**. We hypothesize that Elevation acts through a causal chain
involving Temperature and Soil Moisture:

1.  **Elevation** determines **Temperature** (higher elevation
    $\rightarrow$ lower temperature).
2.  **Temperature** influences **Soil Moisture** (lower temperature
    $\rightarrow$ lower evaporation $\rightarrow$ higher moisture).
3.  **Abundance** is driven by **Moisture**, **Temperature**, and
    potentially a direct effect of **Elevation** (e.g., due to UV
    radiation or partial pressure of gases).

### 1. Simulate Data

We simulate $N = 200$ plots along an elevation gradient.

``` r
library(because)

set.seed(42)
N <- 200

# 1. Elevation (Exogenous variable)
# Ranges roughly from 500m to 1500m
Elevation <- rnorm(N, mean = 1000, sd = 200)

# 2. Temperature (Mediator 1)
# Decreases with Elevation (Lapse rate approx effect)
# Coef: -0.01 implies 100m climb -> -1 degree C
Temp <- 25 - 0.01 * Elevation + rnorm(N, sd = 2)

# 3. Moisture (Mediator 2)
# Driven by Temperature. Cooler -> Moister.
# We model a negative relationship with Temp.
Moisture <- 20 - 2 * Temp + rnorm(N, sd = 5)

# 4. Plant Abundance (Outcome)
# - Positive effect of Moisture (+0.5)
# - Positive effect of Temperature (+1.5)
# - Direct negative effect of Elevation (-0.005) due to harsh conditions
Abundance <- 10 + 0.5 * Moisture + 1.5 * Temp - 0.005 * Elevation + rnorm(N, sd = 10)

eco_data <- data.frame(Elevation, Temp, Moisture, Abundance)
head(eco_data)
```

### 2. Standardize Data

**Important:** For mediation analysis, it is highly recommended to
**standardize** your continuous variables (mean = 0, sd = 1) before
fitting.

Standardization ensures that: 1. **Coefficients are comparable**: All
effects are expressed in standard deviation units (standardized
effects). 2. **Scale Invariance**: The calculation of Indirect Effects
(product of coefficients) is more interpretable relative to the Total
Effect. 3. **Convergence**: MCMC sampling often behaves better with
standardized scales.

``` r
# Standardize all variables
eco_data_std <- scale(eco_data)
head(eco_data_std)
```

### 3. Fit the Structural Equation Model

We define the structural equations reflecting our causal DAG. Notice the
chain: `Elevation -> Temp -> Moisture -> Abundance`.

``` r
# Define the structural equations
eco_eqs <- list(
  Temp ~ Elevation,
  Moisture ~ Temp,
  Abundance ~ Moisture + Temp + Elevation
)

# Fit the model
# We use a short chain for demonstration purposes. Use more iterations for real analysis.
fit <- because(
  equations = eco_eqs,
  data = eco_data_std,
  n.iter = 2000
)
summary(fit)
```

We can also plot the fitted causal model with its standardised paths:

``` r
plot_dag(fit)
```

### 4. Perform Mediation Analysis

We want to understand the Total Effect of **Elevation** on
**Abundance**, and decompose it into its direct and indirect components.

``` r
# Run Mediation Analysis for Elevation -> Abundance
med_results <- because_mediation(fit, exposure = "Elevation", outcome = "Abundance")
```

#### Inspect the Summary

``` r
med_results$summary
```

**Interpretation (Standardized Units):** \* **Total Effect**: The net
impact of Elevation on Abundance in SD units. \* **Direct Effect**: The
path `Elevation -> Abundance` (independent of mediators). \* **Total
Indirect Effect**: The sum of all mediated paths. In our simulation,
Elevation lowers Temp, which raises Moisture, which increases Abundance.

#### Inspect Individual Paths

The `because_mediation` function automatically traces all valid paths
from exposure to outcome in the DAG.

``` r
med_results$paths
```

We expect to see three distinct paths:

1.  **Direct**: `Elevation -> Abundance`
2.  **Short Indirect**: `Elevation -> Temp -> Abundance`
    - Elevation lowers Temp (negative correlation); Temp increases
      Abundance (positive correlation).
    - The product of these effects is **Negative**.
3.  **Long Indirect (Chain)**:
    `Elevation -> Temp -> Moisture -> Abundance`
    - Elevation lowers Temp (negative).
    - Lower Temp raises Moisture (negative relationship $\rightarrow$
      positive change in moisture).
    - Higher Moisture raises Abundance (positive).
    - The chain involves two negative links and one positive link,
      resulting in a **Positive** indirect effect.

The function handles this decomposition automatically, providing
credibility intervals for each specific mechanism.

### Technical Note

The function calculates the indirect effect as the product of
coefficients along the path. For example, for the long chain:
$$\text{Indirect}_{\text{chain}} = \beta_{Elev\rightarrow Temp} \times \beta_{Temp\rightarrow Moist} \times \beta_{Moist\rightarrow Abund}$$
This product-of-coefficients approach is exact for **linear (Gaussian,
identity-link) models** (MacKinnon, 2008; Pearl, 2001).

> **Important warning for non-linear models**: When any equation in the
> mediation path uses a non-linear link function (log for Poisson, logit
> for Bernoulli/Binomial, etc.), the product of link-scale coefficients
> is **not** the natural indirect effect on the outcome scale. The
> distortion can be substantial and can even reverse the direction of
> the indirect effect. For non-linear mediation, use
> counterfactual-based causal mediation methods that operate on the
> response scale, such as those described in VanderWeele (2015) and Imai
> et al. (2010). When working with non-linear outcomes, treat the
> [`because_mediation()`](https://because-pkg.github.io/because/reference/because_mediation.md)
> results as approximate and interpret them cautiously.

### References

Imai, K., Keele, L., & Tingley, D. (2010). A general approach to causal
mediation analysis. *Psychological Methods*, 15(4), 309–334.
<https://doi.org/10.1037/a0020761>

MacKinnon, D. P. (2008). *Introduction to Statistical Mediation
Analysis*. Erlbaum.

Pearl, J. (2001). Direct and indirect effects. *Proceedings of the 17th
Conference on Uncertainty in Artificial Intelligence*, pp. 411–420.
Morgan Kaufmann.

VanderWeele, T. J. (2015). *Explanation in Causal Inference: Methods for
Mediation and Interaction*. Oxford University Press.
