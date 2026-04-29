# Multiscale Causal Models

## Introduction

Many ecological and evolutionary datasets are inherently **multiscale**:
processes operate at different resolutions, such as annual climate
driving individual-level physiology. In such settings, variables
measured at different scales have fundamentally different sample sizes
and thus different statistical power. A climate covariate averaged over
10 years provides only 10 degrees of freedom for testing its effects,
regardless of how many individual observations exist within those years.

The `because` package explicitly handles this structure through its
multiscale engine. When data are supplied as a **named list of data
frames**, one per scale, the engine:

1.  Correctly maps each variable to its natural scale.

2.  Uses the **appropriate degrees of freedom** for each d-separation
    test: year-scale tests use $N_{years}$, individual-scale tests use
    $N_{individuals}$, and cross-scale tests are fitted at the lower
    scale using whatever random effects the researcher specified in the
    structural equations.

3.  Skips tests that are trivially satisfied by the design of the
    multiscale structure (e.g., orthogonal branches that share no
    observations).

## Example: Climate, Resources, and Reproductive Success

Consider a long-term monitoring study of a wild mammal population. Over
10 years, ecologists record:

- **Year-scale** variables (1 row per year): mean winter temperature
  (`Temp`) and annual primary productivity (`NDVI`).
- **Individual-scale** variables (multiple rows per individual per
  year): body mass (`Mass`) and number of offspring (`Offspring`).

The causal hypothesis is:

- Winter temperature controls primary productivity (`Temp → NDVI`).
- Productivity and temperature together drive individual body mass
  (`NDVI → Mass`, `Temp → Mass`).
- Body mass determines reproductive success (`Mass → Offspring`).

This is a **cross-scale** causal chain: climate acts at the year scale,
but its ultimate effect on offspring is mediated through individual body
condition.

``` r
library(because)

dag <- list(
  NDVI       ~ Temp,
  Mass       ~ NDVI + Temp,
  Offspring  ~ Mass
)

plot_dag(dag)
```

## Data Preparation

To use multiscale modeling, your data must be structured as a **named
list of data frames**, where each list element corresponds to a
different scale in your study.

### The “One Home” Rule

Every variable in your model (predictors and responses) must live in
**exactly one** data frame (the one corresponding to the scale where it
was measured).

- **Year-scale variables** (like `Temp`) should only be in the `year`
  data frame.
- **Individual-scale variables** (like `Mass`) should only be in the
  `individual` data frame.

Do **not** repeat year-level variables in the individual-level data
frame. `because` will automatically handle the cross-level mapping using
your `link_vars`.

### Auto-detection vs. Manual Scales

If you omit the `levels` argument, `because` will try to auto-detect
which variable belongs to which scale. This works well for simple
models, but providing the `levels` list manually is recommended
because: 1. It resolves any ambiguity if multiple data frames share
column names. 2. It allows you to explicitly **exclude** nuisance
variables that might be in your data frames but are not part of the
causal model.

## Simulate Multiscale Example Data

``` r
set.seed(42)

n_years <- 10
n_ind   <- 30   # 10 individuals followed over time

# --- Year-scale data ---
year_df <- data.frame(
  Year = 1:n_years,
  Temp = rnorm(n_years, mean = 0, sd = 1),
  NDVI = NA_real_
)
year_df$NDVI <- 0.7 * year_df$Temp + rnorm(n_years, 0, 0.5)

# --- Individual-scale (Repeated Measures) data ---
ind_df <- expand.grid(
  ID   = 1:n_ind,
  Year = 1:n_years
)

# For simulation, each individual has a unique "personality" (random effect)
# that stays constant across years.
ind_personality <- rnorm(n_ind, 0, 0.5)
names(ind_personality) <- 1:n_ind

# Look up year-scale and individual-scale components
temp_lookup <- year_df$Temp[ind_df$Year]
ndvi_lookup <- year_df$NDVI[ind_df$Year]
pers_lookup <- ind_personality[as.character(ind_df$ID)]

# Mass depends on year-scale NDVI/Temp AND the individual's personality
ind_df$Mass <- (
  1.2 * ndvi_lookup + 
  0.4 * temp_lookup + 
  pers_lookup +        # <-- Individual-scale consistency
  rnorm(nrow(ind_df), 0, 0.4)
)

# Offspring depends on Mass
ind_df$Offspring <- rpois(nrow(ind_df), lambda = exp(0.5 + 0.6 * ind_df$Mass))

# Combine into a multiscale list
data_list <- list(
  year       = year_df,
  individual = ind_df
)
```

## Supplying Multiscale Data to `because`

Pass the data as a **named list of data frames**. The names define the
scale labels that `because` uses to assign variables.

``` r
data_list <- list(
  year       = year_df,
  individual = ind_df
)
```

We specify the structural equations as usual:

``` r
equations <- list(
  NDVI      ~ Temp,
  Mass      ~ NDVI + Temp,
  Offspring ~ Mass
)
```

and finally we can fit the model with
[`because()`](https://because-pkg.github.io/because/reference/because.md),
specifying the `multiscale` structure and the variable linking the
datasets with the `link_vars` argument.

``` r
fit <- because(
  equations  = equations,
  data       = data_list,
  family     = c(Offspring = "poisson"),
  multiscale = "year > individual",
  link_vars  = c(year = "Year"),
  n.iter     = 15000,
  n.burnin   = 3000
)

summary(fit)
```

## Understanding the Multiscale Structure

### The `multiscale` argument

`multiscale = "year > individual"` tells `because` that the **year**
scale is coarser and the **individual** scale is nested within it. The
`>` operator means “is a parent of”. (Note: the argument `hierarchy` is
also supported as a synonym).

### The `link_vars` argument

`link_vars = c(year = "Year")` specifies that the column `Year` in the
**individual-scale** data frame is the index that maps each individual
record back to its corresponding row in the **year-scale** data frame.

### Automatic scale assignment

`because` inspects which data frame each variable belongs to:

| Variable    | Scale      | N (df for d-sep) |
|:------------|:-----------|:-----------------|
| `Temp`      | year       | 10               |
| `NDVI`      | year       | 10               |
| `Mass`      | individual | 300              |
| `Offspring` | individual | 300              |

## Random Effects in `because`

While **multiscale** modeling defines the resolution of your variables,
you can still use standard random effects to account for unmeasured
group-level variation. `because` uses the same `(1|grouping_variable)`
notation as `lme4`’s `lmer()` / `glmer()` to specify **random
intercepts**.

A random intercept `(1|Group)` tells `because` that observations sharing
the same value of `Group` are not independent — each group gets its own
intercept drawn from a common Normal distribution.

> \[!NOTE\] Random intercepts specified by the user are propagated into
> the d-separation sub-tests for all equations they apply to. This
> ensures the independence tests reflect exactly the same random
> structure as the structural model, maintaining internal consistency.

### Syntax option A: inline in the equation

``` r
equations_re <- list(
  NDVI      ~ Temp,  
  Mass ~ NDVI + Temp + (1|Year) + (1|ID),    # Year random intercept for Mass
  Offspring ~ Mass  + (1|Year) + (1|ID)  # Year AND individual random intercepts
)
```

### Syntax option B: global `random` argument

``` r
fit_re2 <- because(
  equations  = equations_re,
  data       = data_list,
  family     = c(Offspring = "poisson"),
  multiscale = "year > individual",
  link_vars  = c(year = "Year"),
  random     = ~(1|Year), # applied to ALL equations
  n.iter     = 15000,
  n.burnin   = 3000
)
```

## D-Separation with Correct Degrees of Freedom

The causal model implies the following basis set of conditional
independence claims:

| Test                                    | Scale                  | Conditioning set |
|:----------------------------------------|:-----------------------|:-----------------|
| `NDVI _\|\|_ Mass \| {Temp}`            | **year** (N=10)        | Temp             |
| `NDVI _\|\|_ Offspring \| {Mass, Temp}` | **individual** (N=300) | Mass, Temp       |
| `Temp _\|\|_ Offspring \| {Mass, NDVI}` | **individual** (N=300) | Mass, NDVI       |

> \[!IMPORTANT\] The first test — `NDVI _||_ Mass | {Temp}` — tests a
> **year-scale** independence claim. Even though 300 individual
> observations exist, `because` correctly uses **N=10** degrees of
> freedom for this test, because both NDVI and Temp are year-scale
> variables. Using individual-level N here would be pseudoreplication.

### Within-scale tests

For tests where **both** the response and the focal predictor are at the
**same scale**, `because` fits the sub-model directly on that scale’s
data, naturally using the correct sample size.

### Cross-scale tests

For tests where the response and focal predictor are at **different
scales**, `because` fits the sub-model at the **finer (individual)
scale** using the individual-level data.

Crucially, the d-sep tests use **exactly the same random effect
structure as the structural equations**. If the researcher did not
specify a year-level random effect in their structural model, neither
will the d-sep tests. This is by design and scientifically coherent: the
model assumes that the predictors already in the model fully account for
variation at that scale.

## Visualising the Results

``` r
# DAG with estimated path coefficients
plot_dag(fit)
```

``` r
# D-separation results with 95% credible intervals
plot_dsep(fit)
```

## Practical Guidelines

### When to use multiscale data in `because`

- Your predictors are measured at a **coarser temporal or spatial
  grain** than your response (e.g., annual climate → daily behaviour).
- You have **repeated measures** on the same individuals across multiple
  occasions.
- You want to partition effects into **between-group** (year-scale) and
  **within-group** (individual-scale) components.

### Minimal setup

At minimum, you need:

1.  A **named list** of data frames, one per scale.
2.  A `multiscale` string (e.g., `"year > individual"`).
3.  A `link_vars` named vector mapping the linking column to its parent
    scale.
