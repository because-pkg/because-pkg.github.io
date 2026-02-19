# Deterministic Nodes: Interactions, tresholds and mathematical transformations

In standard Bayesian regression, variables are typically modeled as
**Stochastic Nodes**:
$$Y_{i} \sim \text{Normal}\left( \alpha + \beta X_{i},\sigma \right)$$
This means $Y$ is *correlated* with $X$, but has its own independent
“noise”.

However, in many causal models, some variables are **logical
consequences** of others. For example: **Interactions**:
`Risk = Exposure * Toxicity` **Thresholds**: `IsAdult = (Age > 2)`
**Transformations**: `LogSize = log(Size)`

If you treat `IsAdult` as a stochastic node, the model might learn the
correlation, but it treats the relationship as “soft”. Crucially, if
`Age` is missing (NA), the model might impute `Age = 5` but
`IsAdult = 0`. This violates the laws of logic.

`because` solves this with **Deterministic Nodes**. These are variables
defined by a hard equation (`<-`), ensuring they always stay perfectly
synchronized with their parents, even during imputation and
counterfactual interventions.

## Setup

``` r
library(because)
set.seed(123)
```

## Example 1: Interactions (`A * B`)

Standard R formulas (`y ~ A * B`) handle interactions automatically on
complete data. However, if values of A or B are missing, and unless you
use the standard procedure of eliminating rows with missing data (with a
cosequent loss of power), **imputation** introduces a trap. If you
follow the common workflow of “impute first, then model”, the imputation
step usually assumes simple linear relationships and ignores the
interaction. This often biases the interaction estimate towards zero.

With `because`, you just write the formula naturally. The model imputes
missing values and estimates the interaction **simultaneously**,
ensuring the interaction structure is preserved in the imputed data.

``` r
N <- 100

# Predictors
Rain <- rnorm(N)
Temp <- rnorm(N)

# Interaction: Growth depends on Rain * Temp
Growth <- 0.5 * (Rain * Temp) + rnorm(N, sd = 0.1)

# remove some random values in Rain and Temp
Rain[sample(1:N, 10)] <- NA
Temp[sample(1:N, 15)] <- NA


data <- data.frame(Growth = Growth, Rain = Rain, Temp = Temp)

# Note the formula: Growth ~ Rain * Temp
fit_int <- because(
    equations = list(
        Growth ~ Rain * Temp
    ),
    data = data,
    n.iter = 1000,
    quiet = TRUE
)
# The model creates a deterministic node for the interaction

# Notice the parameter name "beta_Growth_Rain_x_Temp"
summary(fit_int)
#>                          Mean    SD Naive SE Time-series SE   2.5%   50% 97.5%
#> alpha_Growth            0.016 0.023    0.001          0.001 -0.029 0.016 0.060
#> beta_Growth_Rain        0.000 0.025    0.002          0.001 -0.047 0.000 0.049
#> beta_Growth_Rain_x_Temp 0.504 0.024    0.002          0.002  0.455 0.504 0.548
#> beta_Growth_Temp        0.007 0.021    0.001          0.001 -0.034 0.005 0.048
#> sigmaGrowth             0.186 0.014    0.001          0.001  0.160 0.186 0.209
#> sigmaRain               0.932 0.070    0.004          0.004  0.811 0.928 1.084
#> sigmaTemp               0.984 0.070    0.004          0.004  0.864 0.980 1.126
#>                          Rhat n.eff
#> alpha_Growth            0.996   278
#> beta_Growth_Rain        1.013   404
#> beta_Growth_Rain_x_Temp 0.998   240
#> beta_Growth_Temp        0.994   342
#> sigmaGrowth             1.004   314
#> sigmaRain               1.003   240
#> sigmaTemp               1.008   240
#> 
#> DIC:
#> Mean deviance:  375 
#> penalty 25.36 
#> Penalized deviance: 400.3
```

## Example 2: Logic and Thresholds (`Age > 2`)

You can model “tipping points” or categorical definitions. For example,
a marmot (**Marmota marmota**) is only “Mature” if it is older than 2
years.

``` r
# Predictors
Age <- runif(N, 0, 10)

# Deterministic threshold: 1 if Age > 2, else 0
IsMature <- as.numeric(Age > 2)

# Outcome: Mating success depends on Maturity
Mating <- 2 * IsMature + rnorm(N, sd = 0.5)

data_logic <- data.frame(Mating = Mating, Age = Age)

# Use I() to wrap the logic
fit_logic <- because(
    equations = list(
        Mating ~ I(Age > 3)
    ),
    data = data_logic,
    n.iter = 1000,
    quiet = TRUE
)

# The model recovers the effect of the threshold
# Parameter: "beta_Mating_Age_gt_2" where "gt" means "greater than"
summary(fit_logic)
#>                       Mean    SD Naive SE Time-series SE  2.5%   50% 97.5%
#> alpha_Mating         0.363 0.143    0.009          0.010 0.091 0.358 0.639
#> beta_Mating_Age_gt_3 1.663 0.172    0.011          0.012 1.344 1.656 2.005
#> sigmaMating          0.707 0.049    0.003          0.003 0.616 0.707 0.808
#>                       Rhat n.eff
#> alpha_Mating         0.994   215
#> beta_Mating_Age_gt_3 0.994   216
#> sigmaMating          0.998   263
#> 
#> DIC:
#> Mean deviance:  211.9 
#> penalty 3.085 
#> Penalized deviance: 215
```

## Example 3: Multi-Class Definitions

What if you have multiple life stages? \* **Newborn**: Age = 0 \*
**Juvenile**: Age 1 \* **Subadult**: Age 2-4 \* **Adult**: Age \>= 5

You can define this using a sum of logical conditions:

``` r
# Define 4-level class: 1=Newborn, 2=Juv, 3=Sub, 4=Adult
# IMPORTANT: Data must be INTEGERS matching your formula (1, 2, 3, 4), not strings.
# 'Real' Ages are often integers (years), but the process is continuous.
Age <- round(runif(N, 0, 10), 0)

AgeClass <- 1 * (Age == 0) +
    2 * (Age > 0 & Age < 2) +
    3 * (Age >= 2 & Age < 5) +
    4 * (Age >= 5)

# Outcome: Social Rank increases with Life Stage
Rank <- 1.5 * AgeClass + rnorm(N, sd = 0.5)

# Introduce Missing Data in Age
# AgeClass acts as a rigid constraint for imputation.

# Since AgeClass is observed, it 'blocks' information from Rank.

# The model imputes Age primarily to satisfy the deterministic condition (e.g. valid age for the class).

# Note: Input ages are integers, but the model imputes continuous values (e.g. 3.5 years).
Age_True <- Age
Age[1:20] <- NA

data_multi <- data.frame(Rank = Rank, Age = Age, AgeClass = AgeClass)

fit_multi <- because(
    equations = list(
        # 1. Link Age to AgeClass
        # This tells the model: "AgeClass" isn't random; it depends on Age phases.
        AgeClass ~ I(
            1 * (Age < 1) +
                2 * (Age >= 1 & Age < 2) +
                3 * (Age >= 2 & Age < 5) +
                4 * (Age >= 5)
        ),

        # 2. Link AgeClass to Rank (Causal: Class -> Rank)
        Rank ~ AgeClass
    ),
    data = data_multi,
    n.iter = 1000,
    monitor = "all",
    quiet = TRUE
)

# 1. Check effect of AgeClass on Rank
summary(fit_multi)
#>                      Mean    SD Naive SE Time-series SE   2.5%    50%  97.5%
#> alpha_Age           4.836 0.273    0.018          0.016  4.252  4.828  5.395
#> alpha_AgeClass      0.008 0.057    0.004          0.006 -0.100  0.008  0.121
#> alpha_Rank          0.126 0.197    0.013          0.025 -0.276  0.130  0.482
#> beta_Rank_AgeClass  1.455 0.058    0.004          0.007  1.348  1.454  1.577
#> sigmaAge            2.636 0.192    0.012          0.012  2.312  2.630  3.046
#> sigmaAgeClass       0.143 0.009    0.001          0.001  0.126  0.143  0.162
#> sigmaRank           0.498 0.035    0.002          0.002  0.435  0.497  0.573
#> tau_e_Age           0.146 0.021    0.001          0.001  0.108  0.145  0.187
#> tau_e_AgeClass     49.470 6.512    0.420          0.453 37.980 48.949 62.789
#> tau_e_Rank          4.099 0.572    0.037          0.037  3.046  4.050  5.285
#>                     Rhat n.eff
#> alpha_Age          0.998   290
#> alpha_AgeClass     1.009    84
#> alpha_Rank         1.016    68
#> beta_Rank_AgeClass 1.014    71
#> sigmaAge           1.003   240
#> sigmaAgeClass      0.994   274
#> sigmaRank          0.999   240
#> tau_e_Age          1.003   240
#> tau_e_AgeClass     0.994   241
#> tau_e_Rank         0.999   240
#> 
#> DIC:
#> Mean deviance:  313.8 
#> penalty 7.904 
#> Penalized deviance: 321.7

# 2. Inspect the imputed ages
# Notice how the model guesses Age based on the deterministic constraint of AgeClass!
imputed_ages <- extract_imputed(fit_multi)

# Comparison Table
# We compare the Original Class (observed), the True Age (removed), and the Imputed Age
comparison_df <- data.frame(
    OriginalClass = AgeClass[1:20],
    TrueAge = Age_True[1:20],
    ImputedAge = imputed_ages$Mean[1:20]
)
comparison_df
#>    OriginalClass TrueAge ImputedAge
#> 1              2       1  1.5231795
#> 2              4       8  7.0381789
#> 3              4       9  7.0546083
#> 4              3       2  3.6570358
#> 5              4       8  6.9505200
#> 6              4       6  7.2326217
#> 7              4      10  7.0225473
#> 8              3       2  3.5285386
#> 9              4       5  6.9574071
#> 10             4       9  7.2205351
#> 11             4       9  7.0025908
#> 12             1       0 -0.1909213
#> 13             4      10  6.9912593
#> 14             4       5  6.9591981
#> 15             3       4  3.6406078
#> 16             3       4  3.6630148
#> 17             2       1  1.5211577
#> 18             3       2  3.6522039
#> 19             3       4  3.6530379
#> 20             3       3  3.5774251
```

## Example 4: Mathematical transformations (`log(A)`)

You can also use mathematical transformations directly in the formula.

``` r
Mass <- runif(N, 1, 100) # Ensure positive for log
Metabolism <- 0.75 * log(Mass) + rnorm(N, sd = 0.1)

data_math <- data.frame(Metabolism = Metabolism, Mass = Mass)

fit_math <- because(
    equations = list(
        Metabolism ~ log(Mass)
    ),
    data = data_math,
    n.iter = 1000,
    quiet = TRUE
)

# Parameter: "beta_Metabolism_log_Mass"
summary(fit_math)
#>                            Mean    SD Naive SE Time-series SE   2.5%    50%
#> alpha_Metabolism         -0.044 0.069    0.004          0.008 -0.157 -0.047
#> beta_Metabolism_log_Mass  0.762 0.019    0.001          0.002  0.726  0.763
#> sigmaMetabolism           0.177 0.012    0.001          0.001  0.154  0.177
#>                          97.5%  Rhat n.eff
#> alpha_Metabolism         0.096 1.000    80
#> beta_Metabolism_log_Mass 0.798 1.000    86
#> sigmaMetabolism          0.202 0.998   270
#> 
#> DIC:
#> Mean deviance:  -126.4 
#> penalty 3.11 
#> Penalized deviance: -123.3
```

## When do I need `I()`?

You might notice we used [`I()`](https://rdrr.io/r/base/AsIs.html) for
logic and thresholds but not for mathematical transformations. The rule
is:

1.  **Standard Functions**: You interpret things like `log(x)`,
    `sqrt(x)`, `exp(x)` directly as transformations. No
    [`I()`](https://rdrr.io/r/base/AsIs.html) needed.
    - `y ~ log(x)` (OK)
2.  **Interactions**: Use `*` or `:`.
    - `y ~ A * B` (OK)
3.  **Complex Operators**: If you use operators like `>`, `<`, `+`, `-`
    inside a term, R gets confused. Use
    [`I()`](https://rdrr.io/r/base/AsIs.html) to wrap them.
    - `y ~ A + B` (Means “A and B are predictors”)
    - `y ~ I(A + B)` (Means “The SUM of A and B is the predictor”)
    - `y ~ I(A > 5)` (Means “True/False is the predictor”)

## Why “Deterministic” Matters (The Imputation Advantage)

You might ask: *“Why not just create a column `IsMature` before loading
the data?”*

If your data is complete (no NAs), that works fine. But if you have
**missing data**, pre-calculation is dangerous.

### Scenario: The “Broken Logic” of Standard Imputation

Imagine you have a missing `Age` value. You want to impute it.

1.  **Standard Approach (e.g., MICE)**:
    - The imputer sees two separate columns: `Age` (NA) and `IsMature`
      (NA).
    - It guesses `Age = 6` based on body size.
    - It independently guesses `IsMature = 0` based on hormone levels.
    - **Result**: A marmot that is 6 years old but “Immature”.
      **Impossible!**
2.  **The `because` Approach**:
    - `IsMature` is not a variable you impute; it is a **formula**.
    - In MCMC Step 1, the model guesses `Age = 1`. It automatically
      calculates `IsMature = 0`.
    - In MCMC Step 2, the model guesses `Age = 6`. It automatically
      calculates `IsMature = 1`.
    - **Result**: Every single sample in your posterior distribution is
      logically consistent. The “physics” of your model are preserved.

## D-Separation with Deterministic Nodes

How does `because` treat interactions and deterministic nodes in
d-separation tests?

### 1. Standard Interactions (`A * B`)

If you use a standard interaction in your formula:

``` r
equations = list(Y ~ A * B)
```

`because` follows standard DAG theory (Shipley/Pearl). The interaction
describes the **functional form**, not the topology.

- The DAG contains edges: `A -> Y` and `B -> Y`.
- There is **no** separate “interaction node” (`A_x_B`).
- D-separation tests involving `Y` will condition on the parents: `A`
  and `B`.
- You do **not** need to condition on the interaction term itself. Once
  `A` and `B` are known, `A*B` provides no new causal information (it is
  collinear).

### 2. Explicit Deterministic Nodes (`Compound ~ I(A * B)`)

If you explicitly define the interaction as a separate variable (a
deterministic node):

``` r
equations = list(
    Compound ~ I(A * B),
    Y ~ Compound
)
```

`because` treats `Compound` as a distinct node in the DAG.

- The DAG contains edges: `A -> Compound -> Y` and `B -> Compound -> Y`.
- This allows you to test for **Mediation**. You can test if `Y` is
  independent of `A` conditional on `Compound`
  ($Y\bot A \mid Compound$).

### Use Case: The “Blocking” Test

This “Explicit Node” approach is powerful for testing hypotheses about
mechanisms. In Example 3 above, we defined `AgeClass` explicitly:

``` r
AgeClass ~ I(... Age ...)
Rank ~ AgeClass
```

Because `AgeClass` is a distinct node, `because` will test if `Age` is
d-separated from `Rank` given `AgeClass`.

- **Hypothesis**: “Age only affects Rank *through* the Life Stage
  classification.”
- **Test**: `I(Rank, Age | AgeClass)`
- If this test passes, it confirms that `AgeClass` effectively blocks
  the information from `Age`, validating your structural assumption.

``` r
fit_multi_dsep <- because(
    equations = list(
        AgeClass ~ I(
            1 * (Age < 1) +
                2 * (Age >= 1 & Age < 2) +
                3 * (Age >= 2 & Age < 5) +
                4 * (Age >= 5)
        ),
        Rank ~ AgeClass
    ),
    data = data_multi,
    n.iter = 1000,
    dsep = TRUE,
    quiet = TRUE
)

summary(fit_multi_dsep)
#> d-separation Tests
#> ==================
#> 
#> Test: AgeClass _||_ Age | {I(1 * (Age < 1),2 * (Age >= 1 & Age < 2),3 * (Age >= 2 & Age < 5),4 * (Age >= 5))} 
#>          Parameter Estimate LowerCI UpperCI Indep     P  Rhat n.eff
#>  beta_AgeClass_Age   -0.001  -0.017    0.02   Yes 0.942 1.053    51
#> 
#> Test: Rank _||_ Age | {AgeClass} 
#>      Parameter Estimate LowerCI UpperCI Indep    P  Rhat n.eff
#>  beta_Rank_Age   -0.026  -0.084   0.027   Yes 0.35 1.059   120
#> 
#> 
#> Legend:
#>   Indep: 'Yes' = Conditionally Independent, 'No' = Dependent (based on 95% CI)
#>   P: Bayesian probability that the posterior distribution overlaps with zero
```
