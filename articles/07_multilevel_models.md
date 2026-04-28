# Multilevel Causal Models

## Introduction

Many ecological and evolutionary datasets are inherently multilevel:
observations of individuals are nested within years, years within sites,
or individuals within populations. In such settings, variables measured
at different levels of the hierarchy have fundamentally different sample
sizes and thus different statistical power. A climate covariate averaged
over 10 years provides only 10 degrees of freedom for testing its
effects, regardless of how many individual observations exist within
those years.

The `because` package explicitly handles this structure. When data are
supplied as a **named list of data frames**, one per hierarchical level,
the engine:

1.  Correctly maps each variable to its natural level.
2.  Uses the **appropriate degrees of freedom** for each d-separation
    test: year-level tests use $N_{years}$, individual-level tests use
    $N_{individuals}$, and cross-level tests are fitted at the lower
    level using whatever random effects the researcher specified in the
    structural equations.
3.  Skips tests that are trivially satisfied by the design of the
    hierarchy (e.g., orthogonal branches that share no observations).

## A Motivating Example: Climate, Resources, and Reproductive Success

Consider a long-term monitoring study of a wild ungulate population.
Over 10 years, ecologists record:

- **Year-level** variables (1 row per year): mean winter temperature
  (`Temp`), annual primary productivity (`NDVI`).
- **Individual-level** variables (multiple rows per individual per
  year): body mass (`Mass`), number of offspring (`Offspring`).

The causal hypothesis is:

- Winter temperature controls primary productivity (`Temp → NDVI`).
- Productivity and temperature together drive individual body mass
  (`NDVI → Mass`, `Temp → Mass`).
- Body mass determines reproductive success (`Mass → Offspring`).

This is a cross-level causal chain: climate acts at the year level, but
its ultimate effect on offspring is mediated through individual body
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

## Data Preparation and Rules

To use the hierarchical modeling features, your data must be structured
as a **named list of data frames**, where each list element corresponds
to a level in your hierarchy.

### The “One Home” Rule

Every variable in your model (predictors and responses) must live in
**exactly one** data frame—the one corresponding to the level where it
was measured.

- **Year-level variables** (like `Temp`) should only be in the `year`
  data frame.
- **Individual-level variables** (like `Mass`) should only be in the
  `individual` data frame.

Do **not** repeat year-level variables in the individual-level data
frame. `because` will automatically handle the cross-level mapping using
your `link_vars`.

### Auto-detection vs. Manual Levels

If you omit the `levels` argument, `because` will try to auto-detect
which variable belongs to which level. This works well for simple
models, but providing the `levels` list manually is recommended
because: 1. It resolves any ambiguity if multiple data frames share
column names. 2. It allows you to explicitly **exclude** nuisance
variables that might be in your data frames but are not part of the
causal model.

## Step 1: Simulate Hierarchical Data

We simulate a 10-year study with 30 individuals per year. Note how we
keep the year-level and individual-level variables strictly separated.

``` r
set.seed(42)

n_years <- 10
n_ind   <- 30

# --- Year-level data ---
year_df <- data.frame(
  Year = 1:n_years,
  Temp = rnorm(n_years, mean = 0, sd = 1),
  NDVI = NA_real_
)

# NDVI depends on Temp at the year level
year_df$NDVI <- 0.7 * year_df$Temp + rnorm(n_years, 0, 0.5)

# --- Individual-level data ---
ind_df <- data.frame(
  ID         = 1:(n_years * n_ind),
  Year       = rep(1:n_years, each = n_ind),    # <-- linking variable
  Mass       = NA_real_,
  Offspring  = NA_integer_
)

# For simulation, we look up year-level values, 
# but we DO NOT store them in ind_df
temp_lookup <- year_df$Temp[ind_df$Year]
ndvi_lookup <- year_df$NDVI[ind_df$Year]

# Mass depends on year-level NDVI and Temp
ind_df$Mass <- (
  1.2 * ndvi_lookup + 
  0.4 * temp_lookup + 
  rnorm(n_years * n_ind, 0, 0.8)
)

# Offspring depends on individual Mass
ind_df$Offspring <- rpois(n_years * n_ind, lambda = exp(0.5 + 0.6 * ind_df$Mass))

# Combine into a hierarchical list
data_list <- list(
  year       = year_df,
  individual = ind_df
)
```

## Supplying Multilevel Data to `because`

Pass the data as a **named list of data frames**. The names define the
level labels that the engine uses to assign variables.

## Random Effects in `because`

`because` uses the same `(1|grouping_variable)` notation as `lme4`’s
`lmer()` / `glmer()` to specify **random intercepts**. This will be
immediately familiar if you have used mixed-effects models in R.

A random intercept `(1|Group)` tells `because` that observations sharing
the same value of `Group` are not independent — each group gets its own
intercept drawn from a common Normal distribution:

$$u_{g} \sim \text{Normal}\left( 0,\sigma_{u}^{2} \right),\quad g = 1,\ldots,G$$

with the precision $\tau_{u} = 1/\sigma_{u}^{2}$ estimated from the
data.

> \[!NOTE\] **Random slopes** (e.g. `(X|Group)`) are **not yet
> implemented**. Only random intercepts `(1|Group)` are currently
> supported. The syntax deliberately mirrors `lme4` so that random
> slopes can be added as a natural extension in a future release without
> any change to the user interface.

### Syntax: inline in the equation

Place `(1|Group)` directly inside any structural equation, just as you
would in `lmer()`:

``` r
library(because)

# Individuals nested within Sites — Site gets a random intercept in Y
equations <- list(
  Y ~ X + (1|Site),
  Z ~ Y
)

fit <- because(equations = equations, data = my_data)
```

The random term is **equation-specific**: here only `Y` gets the Site
random intercept. `Z ~ Y` is fitted without it (assuming no residual
site-level variation in `Z` beyond what is mediated through `Y`).

You can specify different grouping structures for different responses:

``` r
equations <- list(
  Mass      ~ NDVI + Temp + (1|Year),    # Year random intercept for Mass
  Offspring ~ Mass  + (1|Year) + (1|ID)  # Year AND individual random intercepts
)
```

### Syntax: global `random` argument

When the same random effect applies to **all** response variables, the
`random` argument is a convenient shorthand:

``` r
fit <- because(
  equations = list(Y ~ X, Z ~ Y),
  data      = my_data,
  random    = ~(1|Site)          # (1|Site) added to every equation
)
```

Internally this is identical to embedding `(1|Site)` in each equation.

### When to include a random effect

| Situation                                                                                               | Include `(1|Group)`? |
|:--------------------------------------------------------------------------------------------------------|:---------------------|
| Observations share a year-level environment and year predictors **do not** fully capture year variation | **Yes**              |
| Year-level predictors (e.g. Temp, NDVI) fully explain year variation                                    | Not necessary        |
| Repeated measures on the same individual across time                                                    | **Yes** — `(1|ID)`   |
| Balanced, fully crossed design with all levels as fixed effects                                         | Not needed           |

> \[!IMPORTANT\] Random effects specified in the structural equations
> are automatically carried forward into the d-separation sub-tests for
> all equations they apply to. This keeps the independence tests
> internally consistent with the structural model — see the [Cross-level
> tests](#cross-level-tests) section below.

### Specifying random effects

If you believe there is residual year-level variation beyond what `Temp`
and `NDVI` explain (unmeasured confounders such as snow depth or disease
pressure), you should include a year-level random effect. There are two
equivalent syntaxes — choose the one that fits your workflow:

**Option A — inline in the equation (recommended for equation-specific
effects):**

Just as in `lme4`’s `lmer()`, you can embed `(1|Group)` terms directly
inside the formula. This is the most flexible approach when different
equations require different grouping structures:

``` r
equations_re <- list(
  NDVI      ~ Temp,
  Mass      ~ NDVI + Temp + (1|Year),    # year random intercept for Mass only
  Offspring ~ Mass  + (1|Year)           # year random intercept for Offspring
)

fit_re <- because(
  equations  = equations_re,
  data       = data_list,
  family     = c(Offspring = "poisson"),
  hierarchy  = "year > individual",
  link_vars  = c(year = "Year"),
  dsep       = TRUE,
  n.iter     = 15000,
  n.burnin   = 3000
)
```

**Option B — global `random` argument (convenient when all equations
share the same grouping):**

``` r
fit_re2 <- because(
  equations  = equations,
  data       = data_list,
  family     = c(Offspring = "poisson"),
  hierarchy  = "year > individual",
  link_vars  = c(year = "Year"),
  random     = ~(1|Year),               # applied to ALL equations
  dsep       = TRUE,
  n.iter     = 15000,
  n.burnin   = 3000
)
```

Both syntaxes are internally parsed into the same JAGS model code. The
inline version gives you equation-level control; the global version is a
shorthand when the grouping structure is uniform.

> \[!NOTE\] Random intercepts specified by the user — whether inline or
> global — are propagated into the d-separation sub-tests for all
> equations they apply to. This ensures the independence tests reflect
> exactly the same random structure as the structural model, maintaining
> internal consistency.

For the worked example below we proceed without the year random effect,
assuming `Temp` and `NDVI` fully capture year-level variation:

``` r
data_list <- list(
  year       = year_df,
  individual = ind_df
)

equations <- list(
  NDVI      ~ Temp,
  Mass      ~ NDVI + Temp,
  Offspring ~ Mass
)

fit <- because(
  equations  = equations,
  data       = data_list,
  family     = c(Offspring = "poisson"),
  hierarchy  = "year > individual",
  link_vars  = c(year = "Year"),
  dsep       = TRUE,
  n.iter     = 15000,
  n.burnin   = 3000
)

summary(fit)
```

## Understanding the Hierarchical Data Structure

### The `hierarchy` argument

`hierarchy = "year > individual"` tells `because` that the **year**
level is higher (coarser) and the **individual** level is nested within
it. The `>` operator means “is a parent of”.

### The `link_vars` argument

`link_vars = c(year = "Year")` specifies that the column `Year` in the
**individual-level** data frame is the index that maps each individual
record back to its corresponding row in the **year-level** data frame.
This is analogous to a foreign key in a relational database.

### Automatic level assignment

`because` inspects which data frame each variable belongs to:

| Variable    | Level      | N (df for d-sep) |
|:------------|:-----------|:-----------------|
| `Temp`      | year       | 10               |
| `NDVI`      | year       | 10               |
| `Mass`      | individual | 300              |
| `Offspring` | individual | 300              |

## D-Separation with Correct Degrees of Freedom

The causal model implies the following basis set of conditional
independence claims (the minimal set of claims whose failure would
falsify the model):

| Test                                    | Level                  | Conditioning set |
|:----------------------------------------|:-----------------------|:-----------------|
| `NDVI _\|\|_ Mass \| {Temp}`            | **year** (N=10)        | Temp             |
| `NDVI _\|\|_ Offspring \| {Mass, Temp}` | **individual** (N=300) | Mass, Temp       |
| `Temp _\|\|_ Offspring \| {Mass, NDVI}` | **individual** (N=300) | Mass, NDVI       |

> \[!IMPORTANT\] The first test — `NDVI _||_ Mass | {Temp}` — tests a
> year-level independence claim. Even though 300 individual observations
> exist, `because` correctly uses **N=10** degrees of freedom for this
> test, because both NDVI and Temp (the conditioning set) are year-level
> variables. Using individual-level N here would be pseudoreplication
> and would dramatically inflate statistical power.

### Within-level tests

For tests where **both** the response and the focal predictor are at the
**same level**, `because` fits the sub-model directly on that level’s
data, naturally using the correct sample size.

### Cross-level tests

For tests where the response and focal predictor are at **different
levels** (e.g., testing whether a year-level predictor is independent of
an individual-level response, given a conditioning set), `because` fits
the sub-model at the **lower (individual) level** using the
individual-level data.

Crucially, the d-sep tests use **exactly the same random effect
structure as the structural equations**. If the researcher did not
specify a year-level random effect in their structural model, neither
will the d-sep tests. This is by design and scientifically coherent:

- The structural equation `Offspring ~ Mass` (without `(1|Year)`)
  embodies the scientific assumption that the year-level predictors
  already in the model — `Temp` and `NDVI` — fully account for
  year-level variation in the response. No residual year-level
  heterogeneity needs to be modelled.
- A d-sep test like `Offspring _||_ NDVI | {Mass, Temp}` conditions on
  `Temp`, a year-level variable. Conditioning on `Temp` statistically
  accounts for the year structure *within that test*, making an
  additional `(1|Year)` term redundant given the researcher’s modelling
  assumptions.
- Automatically inserting year random effects into the d-sep tests when
  the researcher deliberately omitted them from the structural model
  would create an **internal inconsistency** — the test would be
  evaluated under stricter assumptions than the model itself.

> \[!IMPORTANT\] If you believe there are **unmeasured year-level
> confounders** beyond `Temp` and `NDVI`, you should include a year
> random effect explicitly in your structural equations (e.g.,
> `random = ~(1|Year)` in
> [`because()`](https://because-pkg.github.io/because/reference/because.md)).
> This will then also be correctly propagated into all relevant d-sep
> tests.

### Orthogonal branches (automatically skipped)

If the hierarchy contains parallel branches (e.g., `site > obs` and
`species > obs`) and the variables in a test come from different
branches with no conditioning variable bridging them, `because` will
skip the test.

This skip is a consequence of the **hierarchical specification**: the
d-separation engine is optimized for linear nesting chains. Tests
between orthogonal branches (which would require crossed random effects)
are currently treated as being outside the structural scope defined by
the linear hierarchy. This ensures the engine remains robust and avoids
attempting cross-level regressions where no shared hierarchical index
exists beyond the observation level.

## Visualising the Results

``` r
# DAG with estimated path coefficients
plot_dag(fit)
```

``` r
# Caterpillar plot of all path coefficients
plot_coef(fit)
```

``` r
# D-separation results with 95% credible intervals
plot_dsep(fit)
```

## Practical Guidelines

### When to use multilevel data in `because`

- Your predictors are measured at a **coarser temporal or spatial
  grain** than your response (e.g., annual climate → daily behaviour).
- You have **repeated measures** on the same individuals across multiple
  occasions (years, seasons).
- You want to partition effects into **between-group** (year-level) and
  **within-group** (individual-level) components.

### Minimal setup

At minimum, you need:

1.  A **named list** of data frames, one per level.
2.  A `hierarchy` string (e.g., `"year > individual"`).
3.  A `link_vars` named vector mapping the linking column to its parent
    level (e.g., `c(year = "Year")`).

If `levels` is omitted, `because` will automatically assign each
variable to the data frame it is found in.

### Variable overlap between levels

If the same column name appears in more than one data frame (e.g., a
site-level mean temperature and an individual-level body temperature
both named `Temp`), use distinct column names to avoid ambiguity.

## References

Shipley, B. (2016). *Cause and Correlation in Biology: A User’s Guide to
Path Analysis, Structural Equations and Causal Inference with R* (2nd
ed.). Cambridge University Press.

Gelman, A., & Hill, J. (2007). *Data Analysis Using Regression and
Multilevel/Hierarchical Models*. Cambridge University Press.

Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
obtaining R² from generalized linear mixed-effects models. *Methods in
Ecology and Evolution*, 4(2), 133–142.
