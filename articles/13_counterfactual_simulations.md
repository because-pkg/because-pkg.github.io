# Counterfactual Simulations with the do() Operator

## Introduction to Counterfactuals

One of the primary goals of causal inference is to answer “what-if”
questions: *What would happen to the population size if we increased the
temperature by 2°C? What if we reduced poaching by 50%?*

In Pearl’s causal framework, these questions are formalized using the
**$`do()`$-calculus**. When you apply an intervention to a variable
(e.g., $`do(X = x)`$), you are essentially forcing that variable to take
a specific value, “severing” it from its natural causes. The effect of
this intervention then ripples downstream through the causal graph.

The `because` package natively supports Bayesian counterfactual
simulations via the
[`do()`](https://because-pkg.github.io/because/reference/do.md)
operator. Because the engine understands the structural equations and
the topology of your DAG, it automatically propagates interventions
correctly while carrying forward the full posterior uncertainty of the
model.

``` r

library(because)
```

## Simulating a Causal System

Let’s set up a classic ecological cascade: **Temperature** drives the
availability of a **Resource**, which in turn determines the
**Abundance** of a species. We will also include a direct physiological
effect of Temperature on Abundance.

``` r

set.seed(42)
N <- 100

# True parameters
# Temperature (mean = 20C, sd = 4C)
raw_Temp <- rnorm(N, mean = 20, sd = 4)

# We scale variables for model fitting (Standard practice)
df <- data.frame(Site = 1:N)
df$Temp <- scale(raw_Temp) 

# Resource is driven by Temperature
df$Resource <- 0.6 * df$Temp + rnorm(N, 0, 1)

# Abundance is driven by both (Poisson distributed)
log_lambda <- 1.2 * df$Resource - 0.5 * df$Temp
df$Abundance <- rpois(N, lambda = exp(log_lambda))
```

We define our equations and fit the structural model:

``` r

eqs <- list(
  Resource ~ Temp,
  Abundance ~ Resource + Temp
)

# Fit the model using JAGS (short chains for vignette speed)
fit <- because(eqs, data = df, 
               family = list(Abundance = "poisson"),
               n.iter = 500, n.burnin = 100, n.chains = 2, 
               quiet = TRUE)
```

## The `do()` Operator

The [`do()`](https://because-pkg.github.io/because/reference/do.md)
function takes a fitted `because` model and a set of targeted
interventions. It returns a `because_counterfactual` object containing
matrices of simulated posterior draws for every variable in the graph.

### 1. Atomic (Hard) Interventions

An atomic intervention forces a variable to a specific, fixed value
across all observations. Because we scaled our data, `Temp = 2` means we
are setting the temperature to exactly +2 standard deviations above the
historical mean for every single site.

``` r

# Force Temperature to +2 SD at all sites
res_atomic <- do(fit, Temp = 2)

# Summarize the counterfactual predictions globally
summary(res_atomic)
#> Counterfactual Simulation Summary
#> ---------------------------------
#> Estimates represent the global expectation (averaged across all observations)
#> under the intervened causal structure.
#> 
#>   Variable  Mean     SD   2.5%   50% 97.5%
#>       Temp 2.000 0.0000 2.0000 2.000 2.000
#>   Resource 1.147 0.2395 0.7416 1.161 1.538
#>  Abundance 2.143 0.7733 1.0383 2.075 3.803
```

Notice that the `SD` for `Temp` is exactly 0—because we intervened and
forced it to be deterministically fixed. However, `Resource` and
`Abundance` maintain their Bayesian uncertainty.

### 2. Shift (Additive) Interventions

Often, we don’t want to force all sites to be identical. Instead, we
want to shift the historical baseline. We can do this using a formula,
where `.` represents the historical/natural values of the variable.

``` r

# Increase historical temperature by 1 standard deviation at every site
res_shift <- do(fit, Temp = ~ . + 1)
summary(res_shift)
#> Counterfactual Simulation Summary
#> ---------------------------------
#> Estimates represent the global expectation (averaged across all observations)
#> under the intervened causal structure.
#> 
#>   Variable  Mean     SD   2.5%    50%  97.5%
#>       Temp 1.000 0.0000 1.0000 1.0000 1.0000
#>   Resource 0.537 0.1841 0.1228 0.5155 0.8798
#>  Abundance 1.832 0.4691 1.1080 1.8100 2.7860
```

## Intervening on the Raw Metric (`raw_scale = TRUE`)

There is a major conceptual trap when applying interventions to z-scored
(scaled) data. As we just saw, adding `+1` to scaled data means adding 1
standard deviation, not 1 real-world unit.

If you assigned your scaled data safely (`df$Temp <- scale(raw_Temp)`),
`because` secretly saved the original raw mean and standard deviation.
By using `raw_scale = TRUE`, you can tell the
[`do()`](https://because-pkg.github.io/because/reference/do.md) operator
to temporarily unscale the data, apply your intervention on the **raw
biological metric** (like +2°C), and then re-scale it back before
running the simulation!

``` r

# Increase the RAW temperature by exactly 2 degrees Celsius
# (The model handles the standard deviation math for you)
res_raw <- do(fit, Temp = ~ . + 2, raw_scale = TRUE)

summary(res_raw)
#> Counterfactual Simulation Summary
#> ---------------------------------
#> Estimates represent the global expectation (averaged across all observations)
#> under the intervened causal structure.
#> 
#>   Variable   Mean     SD  2.5%   50%  97.5%
#>       Temp 22.130 0.0000 22.13 22.13 22.130
#>   Resource 21.073 0.5682 20.06 21.10 22.252
#>  Abundance  1.661 0.4414  1.06  1.57  2.532
```

> **Important Note:** Base R’s `data.frame(Temp = scale(raw_temp))`
> silently strips scaling attributes. To use `raw_scale = TRUE`, you
> must add scaled columns to existing dataframes using the `$` operator:
> `df$Temp <- scale(raw_temp)`.

### 3. Percentage Shifts

You can also apply multiplicative shifts easily by passing a character
string ending in `%`.

**However, remember our warning about scaled data!** If a variable has a
mean of 0, multiplying it by 1.10 just stretches the variance. To safely
apply a percentage shift, you should use the `raw_scale = TRUE`
argument. This tells the engine to unscale the data to its natural
positive values, apply the 10% increase, and automatically rescale it
back before simulation.

``` r

# Increase the RAW historical temperature by 10% 
# (e.g., if it was 20C, it becomes 22C)
res_perc <- do(fit, Temp = "+10%", raw_scale = TRUE)
summary(res_perc)
#> Counterfactual Simulation Summary
#> ---------------------------------
#> Estimates represent the global expectation (averaged across all observations)
#> under the intervened causal structure.
#> 
#>   Variable   Mean     SD   2.5%    50%  97.5%
#>       Temp 22.143 0.0000 22.143 22.143 22.143
#>   Resource 21.024 0.5649 19.847 21.147 22.001
#>  Abundance  1.619 0.4323  1.089  1.535  2.387
```

*(Note: Under the hood, the `%` string syntax is a convenient shortcut
that intercepts the unscaled matrix and multiplies it by
`1 + (val/100)`).*

### 4. Stochastic Interventions

Sometimes policies are not perfectly exact. You can introduce noise into
your intervention using standard R random number generators. Use `n` in
the formula to represent the number of data points.

``` r

# Shift temperature by +1 SD, but add observation noise to the intervention
res_stoch <- do(fit, Temp = ~ rnorm(n, mean = . + 1, sd = 0.2))
summary(res_stoch)
#> Counterfactual Simulation Summary
#> ---------------------------------
#> Estimates represent the global expectation (averaged across all observations)
#> under the intervened causal structure.
#> 
#>   Variable   Mean      SD   2.5%    50%  97.5%
#>       Temp 1.0031 0.02076 0.9617 1.0027 1.0409
#>   Resource 0.5548 0.16621 0.2813 0.5637 0.8201
#>  Abundance 1.8420 0.52583 1.1093 1.7550 3.0337
```

## Extracting Site-Specific Counterfactuals

The [`summary()`](https://rdrr.io/r/base/summary.html) function is great
for looking at the **global expected value** (averaged across all
sites). But often, you want to map or analyze how specific sites
responded to the intervention.

The [`do()`](https://because-pkg.github.io/because/reference/do.md)
operator returns the full `[ndraws x N_obs]` matrices.

``` r

# Extract the posterior matrix for Abundance under the raw shift
abund_matrix <- res_raw$Abundance

# Look at the dimensions: [posterior draws, sites]
dim(abund_matrix)
#> [1]  80 100

# Calculate the mean counterfactual abundance for Site 5
mean(abund_matrix[, 5])
#> [1] 1.7

# Calculate the 95% Credible Interval for Site 5
quantile(abund_matrix[, 5], probs = c(0.025, 0.975))
#>  2.5% 97.5% 
#> 0.000 9.125
```

This flexibility allows you to easily compute site-specific treatment
effects (the difference between the historical values and the
counterfactual values) for precision ecology!
