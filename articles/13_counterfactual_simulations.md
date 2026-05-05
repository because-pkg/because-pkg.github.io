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
simulations via the `do()` operator. Because the engine understands the
structural equations and the topology of your DAG, it automatically
propagates interventions correctly while carrying forward the full
posterior uncertainty of the model.

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

The `do()` function takes a fitted `because` model and a set of targeted
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
#>   Resource 1.138 0.2197 0.7211 1.119 1.578
#>  Abundance 2.139 0.9129 1.1583 1.910 4.306
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
#>   Resource 0.5327 0.1566 0.2635 0.5329 0.9146
#>  Abundance 1.8322 0.6464 1.0675 1.7000 3.4912
```

## Intervening on the Raw Metric (`raw_scale = TRUE`)

There is a major conceptual trap when applying interventions to z-scored
(scaled) data. As we just saw, adding `+1` to scaled data means adding 1
standard deviation, not 1 real-world unit.

If you assigned your scaled data safely (`df$Temp <- scale(raw_Temp)`),
`because` secretly saved the original raw mean and standard deviation.
By using `raw_scale = TRUE`, you can tell the `do()` operator to
temporarily unscale the data, apply your intervention on the **raw
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
#>   Variable   Mean     SD   2.5%    50%  97.5%
#>       Temp 22.130 0.0000 22.130 22.130 22.130
#>   Resource 21.047 0.5629 19.754 21.076 21.996
#>  Abundance  1.655 0.5526  1.009  1.525  2.601
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
#>   Variable   Mean     SD  2.5%   50%  97.5%
#>       Temp 22.143 0.0000 22.14 22.14 22.143
#>   Resource 21.001 0.4816 20.09 21.06 21.886
#>  Abundance  1.601 0.4172  1.01  1.52  2.696
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
#>       Temp 1.0030 0.02077 0.9617 1.0019 1.0409
#>   Resource 0.5491 0.16099 0.2024 0.5413 0.8858
#>  Abundance 1.8321 0.75153 0.9690 1.6800 4.0965
```

## Extracting Site-Specific Counterfactuals

The [`summary()`](https://rdrr.io/r/base/summary.html) function is great
for looking at the **global expected value** (averaged across all
sites). But often, you want to map or analyze how specific sites
responded to the intervention.

The `do()` operator returns the full `[ndraws x N_obs]` matrices.

``` r

# Extract the posterior matrix for Abundance under the raw shift
abund_matrix <- res_raw$Abundance

# Look at the dimensions: [posterior draws, sites]
dim(abund_matrix)
#> [1]  80 100

# Calculate the mean counterfactual abundance for Site 5
mean(abund_matrix[, 5])
#> [1] 1.7875

# Calculate the 95% Credible Interval for Site 5
quantile(abund_matrix[, 5], probs = c(0.025, 0.975))
#>  2.5% 97.5% 
#> 0.000 7.075
```

This flexibility allows you to easily compute site-specific treatment
effects (the difference between the historical values and the
counterfactual values) for precision ecology!
