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
#> alpha_Growth            0.015 0.011    0.001          0.001 -0.007 0.014 0.039
#> beta_Growth_Rain        0.001 0.014    0.001          0.001 -0.023 0.002 0.025
#> beta_Growth_Rain_x_Temp 0.508 0.013    0.001          0.001  0.480 0.509 0.531
#> beta_Growth_Temp        0.009 0.010    0.001          0.001 -0.011 0.009 0.030
#> sigmaGrowth             0.096 0.008    0.001          0.001  0.083 0.096 0.113
#> sigmaRain               0.945 0.066    0.004          0.004  0.820 0.944 1.084
#> sigmaTemp               0.977 0.073    0.005          0.004  0.853 0.974 1.151
#> sigma_e_Growth          0.096 0.008    0.001          0.001  0.083 0.096 0.113
#> sigma_e_Rain            0.945 0.066    0.004          0.004  0.820 0.944 1.084
#> sigma_e_Temp            0.977 0.073    0.005          0.004  0.853 0.974 1.151
#>                          Rhat n.eff
#> alpha_Growth            1.001   237
#> beta_Growth_Rain        1.006   240
#> beta_Growth_Rain_x_Temp 1.006   187
#> beta_Growth_Temp        0.997   467
#> sigmaGrowth             0.995   246
#> sigmaRain               0.998   240
#> sigmaTemp               1.005   359
#> sigma_e_Growth          0.995   246
#> sigma_e_Rain            0.998   240
#> sigma_e_Temp            1.005   359
#> 
#> DIC:
#> Mean deviance:  300.1 
#> penalty 29.02 
#> Penalized deviance: 329.2
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
#> alpha_Mating         0.356 0.129    0.008          0.009 0.111 0.352 0.618
#> beta_Mating_Age_gt_3 1.667 0.157    0.010          0.010 1.401 1.656 1.975
#> sigmaMating          0.710 0.053    0.003          0.003 0.619 0.706 0.813
#> sigma_e_Mating       0.710 0.053    0.003          0.003 0.619 0.706 0.813
#>                       Rhat n.eff
#> alpha_Mating         1.005   209
#> beta_Mating_Age_gt_3 1.004   240
#> sigmaMating          0.998   272
#> sigma_e_Mating       0.998   272
#> 
#> DIC:
#> Mean deviance:  212.1 
#> penalty 3.252 
#> Penalized deviance: 215.3
```

## Example 3: Multi-Class Definitions

What if you have multiple life stages? \* **Newborn**: Age = 0 \*
**Juvenile**: Age 1 \* **Subadult**: Age 2-4 \* **Adult**: Age \>= 5

You can define this using a sum of logical conditions:

``` r
# Define 4-level class: 1=Newborn, 2=Juv, 3=Sub, 4=Adult
# IMPORTANT: The class must be derived from Age so the model understands the structure.
Age <- round(runif(N, 0, 10), 0)

AgeClass <- 1 * (Age == 0) +
    2 * (Age > 0 & Age < 2) +
    3 * (Age >= 2 & Age < 5) +
    4 * (Age >= 5)

# Outcome: Social Rank increases with Life Stage
Rank <- 1.5 * AgeClass + rnorm(N, sd = 0.5)

data_multi <- data.frame(Rank = Rank, Age = Age, AgeClass = AgeClass)

# because() recovers the slope (~1.5) and confirms the causal path Age -> AgeClass -> Rank
fit_multi <- because(
    equations = list(
        # 1. Link AgeClass to Rank (Causal: Class -> Rank)
        Rank ~ AgeClass
    ),
    data = data_multi,
    n.iter = 1000,
    quiet = TRUE
)

# Check effect of AgeClass on Rank
summary(fit_multi)
#>                     Mean    SD Naive SE Time-series SE   2.5%   50% 97.5%  Rhat
#> alpha_Rank         0.126 0.200    0.013          0.035 -0.257 0.136 0.497 0.997
#> beta_Rank_AgeClass 1.457 0.057    0.004          0.010  1.354 1.456 1.570 0.998
#> sigmaRank          0.473 0.037    0.002          0.002  0.413 0.470 0.566 0.996
#> sigma_e_Rank       0.473 0.037    0.002          0.002  0.413 0.470 0.566 0.996
#>                    n.eff
#> alpha_Rank            38
#> beta_Rank_AgeClass    43
#> sigmaRank            226
#> sigma_e_Rank         226
#> 
#> DIC:
#> Mean deviance:  134.6 
#> penalty 2.941 
#> Penalized deviance: 137.5
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
#> alpha_Metabolism         -0.050 0.041    0.003          0.005 -0.126 -0.049
#> beta_Metabolism_log_Mass  0.764 0.011    0.001          0.001  0.744  0.764
#> sigmaMetabolism           0.105 0.008    0.000          0.000  0.093  0.104
#> sigma_e_Metabolism        0.105 0.008    0.000          0.000  0.093  0.104
#>                          97.5%  Rhat n.eff
#> alpha_Metabolism         0.027 1.018    75
#> beta_Metabolism_log_Mass 0.786 1.022    83
#> sigmaMetabolism          0.121 0.998   240
#> sigma_e_Metabolism       0.121 0.998   240
#> 
#> DIC:
#> Mean deviance:  -166.7 
#> penalty 3.077 
#> Penalized deviance: -163.6
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
#>   (Safety) Dropping cyclic d-sep test: AgeClass ~ Age | I(1 * (Age < 1) + 2 * (Age >= 1 & Age < 2) + 3 * (Age >= 2 & Age < 5) + 4 * (Age >= 5)) (conditioned on own component)

summary(fit_multi_dsep)
#> d-separation Tests
#> ==================
#> 
#> Test: Rank _||_ Age | {AgeClass} 
#>      Parameter Estimate LowerCI UpperCI  Rhat n.eff
#>  beta_Rank_Age   -0.005  -0.072   0.067 1.093    49
```
