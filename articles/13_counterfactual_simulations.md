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
#>   Variable  Mean     SD  2.5%   50% 97.5%
#>       Temp 2.000 0.0000 2.000 2.000 2.000
#>   Resource 1.141 0.2516 0.680 1.105 1.628
#>  Abundance 2.127 0.7548 1.189 1.880 4.054
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
#>   Variable   Mean     SD   2.5%    50%  97.5%
#>       Temp 1.0000 0.0000 1.0000 1.0000 1.0000
#>   Resource 0.5377 0.1707 0.1911 0.5421 0.8787
#>  Abundance 1.8226 0.4966 1.0595 1.7200 3.1122
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
#>   Variable   Mean     SD   2.5%   50% 97.5%
#>       Temp 22.130 0.0000 22.130 22.13 22.13
#>   Resource 21.074 0.6114 19.999 21.08 22.33
#>  Abundance  1.641 0.4063  1.018  1.58  2.74
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
#>   Variable   Mean     SD   2.5%   50%  97.5%
#>       Temp 22.143 0.0000 22.143 22.14 22.143
#>   Resource 21.046 0.6516 19.799 21.08 22.018
#>  Abundance  1.614 0.4326  1.048  1.58  2.637
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
#>       Temp 0.9993 0.02104 0.9523 1.0038 1.0347
#>   Resource 0.5344 0.17032 0.2598 0.5363 0.7735
#>  Abundance 1.7719 0.50792 1.0398 1.6600 2.7905
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
#> [1] 1.875

# Calculate the 95% Credible Interval for Site 5
quantile(abund_matrix[, 5], probs = c(0.025, 0.975))
#>  2.5% 97.5% 
#>  0.00  8.05
```

This flexibility allows you to easily compute site-specific treatment
effects (the difference between the historical values and the
counterfactual values) for precision ecology!
