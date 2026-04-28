# Latent Variables and the MAG Approach

## Introduction

In causal inference, a **latent variable** is a variable that influences
the system but was never measured. Unmeasured common causes are
particularly problematic because they create spurious correlations
between the variables they affect, which can be mistaken for direct
causal effects.

The standard DAG framework assumes all common causes are observed. When
there are unmeasured common causes, the appropriate graphical model is a
**Maximal Ancestral Graph (MAG)**, which represents unobserved common
causes with **bidirected edges** (e.g.,
$\left. X\leftrightarrow Y \right.$).

`because` supports latent variables through the `latent` argument, and
adapts the d-separation basis set accordingly using the MAG framework of
Shipley & Douda (2021).

------------------------------------------------------------------------

### Background: DAGs, MAGs, and Bidirected Edges

In a standard Directed Acyclic Graph (DAG), every causal relationship is
represented by a directed edge ($\left. X\rightarrow Y \right.$). All
variables — observed and unobserved — appear as nodes.

In a **Maximal Ancestral Graph (MAG)**, latent variables are
**marginalized out**, leaving only the observed variables. The effect of
a latent common cause $L$ of two observed variables $X$ and $Y$ is
represented by a **bidirected edge**
$\left. X\leftrightarrow Y \right.$, indicating that the residuals of
$X$ and $Y$ are correlated.

This matters for d-separation: conditioning on different variables has
different consequences in a MAG compared to a DAG, and the basis set of
conditional independence claims must be derived from the MAG, not the
original DAG.

------------------------------------------------------------------------

### The `latent` Argument

Declare latent variables by passing their names to the `latent` argument
of
[`because()`](https://because-pkg.github.io/because/reference/because.md).
The model must still include **equations for the latent variable** (so
its effects on observed variables are estimated), but the latent
variable itself will not be present in the data. `because` will treat it
as an unobserved node and compute the MAG-based d-separation basis set
automatically.

------------------------------------------------------------------------

### Example: Unmeasured Individual Quality

#### Study Design

We monitor 150 individual birds over a breeding season. We measure:

- **Food_Availability** (X₁): Territory food resources (observed).
- **Clutch_Size** (Y₁): Number of eggs laid (observed).
- **Fledgling_Success** (Y₂): Number of young fledged (observed).

We hypothesize that an unmeasured **Individual_Quality** (L) — a
combination of genetic quality, immune status, and experience — affects
both `Clutch_Size` and `Fledgling_Success` simultaneously. Without
accounting for this latent variable, the correlation between
`Clutch_Size` and `Fledgling_Success` might be mistakenly attributed to
a direct causal path between them.

#### The Causal Structure

The underlying DAG (with the latent variable) is:

    Food_Availability → Clutch_Size
    Food_Availability → Fledgling_Success
    Individual_Quality → Clutch_Size
    Individual_Quality → Fledgling_Success

When we marginalise over `Individual_Quality`, the MAG contains a
bidirected edge:

    Food_Availability → Clutch_Size
    Food_Availability → Fledgling_Success
    Clutch_Size ↔ Fledgling_Success   (residual correlation from latent quality)

#### Simulate Data

``` r
library(because)

set.seed(42)
N <- 150

# Observed predictor
Food_Availability <- rnorm(N, 0, 1)

# Latent variable (unobserved — we simulate it for data generation only)
Individual_Quality <- rnorm(N, 0, 1)

# Outcomes driven by both observed and latent causes
Clutch_Size <- round(
  4 + 0.8 * Food_Availability + 1.2 * Individual_Quality + rnorm(N, 0, 0.5)
)
Clutch_Size <- pmax(Clutch_Size, 1)  # Counts must be positive

Fledgling_Success <- round(
  2 + 0.6 * Food_Availability + 0.9 * Individual_Quality + rnorm(N, 0, 0.5)
)
Fledgling_Success <- pmax(Fledgling_Success, 0)

# Data frame — Individual_Quality is NOT included (it's unobserved)
bird_data <- data.frame(
  Food_Availability  = scale(Food_Availability)[, 1],
  Clutch_Size        = Clutch_Size,
  Fledgling_Success  = Fledgling_Success
)
```

#### Fit the Model with a Latent Variable

``` r
# Structural equations: the latent variable appears on the RHS of both equations.
# because() detects that Individual_Quality is not in the data and treats it as latent.
equations_latent <- list(
  Clutch_Size       ~ Food_Availability + Individual_Quality,
  Fledgling_Success ~ Food_Availability + Individual_Quality
)

fit_latent <- because(
  equations = equations_latent,
  data      = bird_data,
  latent    = "Individual_Quality",  # Declare the latent variable
  dsep      = TRUE,                  # Uses MAG-based d-sep basis set
  n.iter    = 20000,
  n.burnin  = 5000
)

summary(fit_latent)
```

The summary now includes estimates for:

- `beta_Clutch_Size_Food_Availability`: effect of food on clutch size.
- `beta_Fledgling_Success_Food_Availability`: effect of food on
  fledgling success.
- `beta_Clutch_Size_Individual_Quality`: estimated loading of the latent
  quality on clutch size.
- `beta_Fledgling_Success_Individual_Quality`: estimated loading on
  fledgling success.
- A **residual correlation** between `Clutch_Size` and
  `Fledgling_Success` that captures the shared latent influence.

#### D-Separation with a MAG

When `latent` is specified, `because` automatically derives the
conditional independence basis set from the **MAG** rather than the DAG.
The MAG-based tests correctly account for the bidirected edge between
`Clutch_Size` and `Fledgling_Success`.

Without the latent variable, the DAG would incorrectly imply that
`Clutch_Size ⊥ Fledgling_Success | {Food_Availability}` — a test that
would **fail** in the data because of the shared latent influence. With
the MAG representation, this test is dropped from the basis set
(bidirected edges imply residual dependence even after conditioning on
observed parents), giving an honest assessment of model fit.

``` r
# Inspect d-separation results (MAG-based basis set)
# The tests shown reflect the independence claims implied by the MAG,
# not the original full DAG.
print(fit_latent$dsep_results)
```

------------------------------------------------------------------------

### Controlling Latent Variable Behaviour

#### `latent_method`

The `latent_method` argument controls how the latent variable is
integrated into the model:

- **`"correlations"`** (default): Introduces a residual correlation
  parameter between the observed descendants of the latent variable.
  Equivalent to marginalising over the latent and estimating a free
  covariance.
- **`"factor"`**: Models the latent variable as an explicit factor,
  estimating loadings and a precision. More interpretable but requires
  identifiability constraints (see `fix_latent`).

``` r
fit_factor <- because(
  equations  = equations_latent,
  data       = bird_data,
  latent     = "Individual_Quality",
  latent_method = "factor",       # Explicit factor model
  fix_latent    = "loading",      # Fix first loading to 1 for identifiability
  standardize_latent = TRUE,      # Standardise the latent factor (SD = 1)
  n.iter     = 20000,
  n.burnin   = 5000
)

summary(fit_factor)
```

#### `fix_latent`

When using `latent_method = "factor"`, the factor scale is not
identified without a constraint. `fix_latent = "loading"` fixes the
first loading (the first equation’s latent effect) to 1, anchoring the
latent scale.

#### `standardize_latent`

Setting `standardize_latent = TRUE` standardises the latent variable to
have unit variance, making the loadings comparable across variables and
easier to interpret as standardised path coefficients.

------------------------------------------------------------------------

### Comparing Models: With and Without the Latent Variable

``` r
# Model 1: No latent variable (misspecified — direct path between outcomes)
fit_no_latent <- because(
  equations = list(
    Clutch_Size       ~ Food_Availability,
    Fledgling_Success ~ Food_Availability + Clutch_Size  # Direct path!
  ),
  data   = bird_data,
  WAIC   = TRUE,
  quiet  = TRUE
)

# Model 2: Latent variable (correct specification)
fit_with_latent <- because(
  equations  = equations_latent,
  data       = bird_data,
  latent     = "Individual_Quality",
  WAIC       = TRUE,
  quiet      = TRUE
)

because_compare(fit_no_latent, fit_with_latent)
```

The latent variable model should have a lower WAIC when the data truly
contain a shared unmeasured cause.

------------------------------------------------------------------------

### Practical Notes

#### When to use a latent variable

- You have **theoretical reasons** to expect a shared unmeasured cause
  (e.g., genetic quality, unobserved environmental heterogeneity).
- A d-separation test **fails** when you know the causal structure is
  correct — a latent common cause may be creating residual dependence.
- A direct path between two outcomes makes biological sense only if you
  include the latent common cause.

#### Identifiability

Latent variable models are only identified when the latent variable has
at least **two observed indicator variables** (i.e., it must influence
at least two observed outcomes). Models with a single indicator are not
identified.

#### Sample size

Estimating latent variables requires more data than purely observed
models. With small samples ($N < 50$), the latent loadings may have very
wide credible intervals.

------------------------------------------------------------------------

### References

Shipley, B., & Douma, J. C. (2021). A new test for d-separation in
maximal ancestral graphs. *Structural Equation Modeling*, 28(2),
167–180.

Pearl, J. (2009). *Causality: Models, Reasoning, and Inference* (2nd
ed.). Cambridge University Press.

Richardson, T., & Spirtes, P. (2002). Ancestral graph Markov models.
*The Annals of Statistics*, 30(4), 962–1030.

von Hardenberg, A., & Gonzalez-Voyer, A. (2025). PhyBaSE: Phylogenetic
Bayesian Structural Equation modelling. *Methods in Ecology and
Evolution*.
