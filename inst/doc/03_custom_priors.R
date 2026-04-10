## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(because)

set.seed(42)
N <- 30
# Ambient temperature in Celsius (Mean 20C)
Temp_Raw <- rnorm(N, mean = 20, sd = 5)

# Center it so '0' represents the mean (20C).
# This ensures the intercept represents growth at average conditions,
# rather than at 0°C (where the bird would be frozen solid!).
Temp_Centered <- Temp_Raw - 20

# True relationship:
# Baseline growth at mean temp is 10 g/day.
# Every degree of warming above mean adds 0.5 g/day.
Growth_g_day <- 0.5 * Temp_Centered + 10 + rnorm(N, sd = 2)

df <- data.frame(Temp_Centered, Growth_g_day)

## ----default_model------------------------------------------------------------
fit_default <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df
)

summary(fit_default)

## ----custom_model-------------------------------------------------------------
# Define our custom priors
my_priors <- list(
    # Strong prior on intercept: Mean 10, Precision 100 (SD = 0.1)
    alphaGrowth_g_day = "dnorm(10, 100)",

    # Informative prior on slope: Mean 0.5, Precision 400 (SD = 0.05)
    beta_Growth_g_day_Temp_Centered = "dnorm(0.5, 400)"
)

fit_custom <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    priors = my_priors
)

summary(fit_custom)

## ----plot_compare, fig.width=9, fig.height=4----------------------------------
# We can use the helper function plot_posterior() to compare models.
# By passing a list of models, they are overlaid on the same plot.
# The 'parameter' argument uses partial matching (regex), so "Growth_g_day"
# will match both the intercept (alphaGrowth_g_day) and the slope (beta_Growth_g_day...).

plot_posterior(
    model = list("Default" = fit_default, "Custom" = fit_custom),
    parameter = "Growth_g_day"
)

## ----mechanistic_constraint---------------------------------------------------
# Simulate data following Kleiber's Law: MR = a * Mass^0.75
# Taking logs: log(MR) = log(a) + 0.75 * log(Mass)
set.seed(42)
Mass <- runif(30, 10, 100)
# True scaling exponent is 0.75
# We add enough noise that the estimated slope might be negative by chance
Log_Mass <- log(Mass)
Log_MR <- 1.5 + 0.75 * Log_Mass + rnorm(30, sd = 2.5)

df_kleiber <- data.frame(Log_Mass, Log_MR)

# Prior: The scaling exponent (slope) must be positive
# Physics dictates that metabolic cost increases with mass.
positive_prior <- list(
    beta_Log_MR_Log_Mass = "dnorm(0, 1) T(0, )"
)

# 1. Fit Default Model (Unconstrained) - Log-Log regression
fit_default_kleiber <- because(
    equations = list(Log_MR ~ Log_Mass),
    data = df_kleiber
)

# 2. Fit Mechanistic Model (Constrained)
fit_mech_kleiber <- because(
    equations = list(Log_MR ~ Log_Mass),
    data = df_kleiber,
    priors = positive_prior
)

# Compare Estimates
summary(fit_default_kleiber)
summary(fit_mech_kleiber)

# Visualize: Unconstrained vs. Truncated
plot_posterior(
    list(Unconstrained = fit_default_kleiber, Mechanistic = fit_mech_kleiber),
    parameter = "beta_Log_MR_Log_Mass",
    density_args = list(Mechanistic = list(from = 0))
)

## ----variance_prior-----------------------------------------------------------
# Simulate small, noisy dataset
set.seed(42)
N_small <- 15
Temp_Small <- rnorm(N_small, 20, 5)
# Actual residual SD = 3 (quite noisy)
Growth_Small <- 0.5 * Temp_Small + rnorm(N_small, sd = 3)
df_small <- data.frame(Growth_g_day = Growth_Small, Temp_Centered = Temp_Small - 20)

# 1. Default Model
# dunif(0, 100) on sigma_e: completely open to any residual SD up to 100.
# With only 15 observations the posterior will be wide.
fit_default_var <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df_small,
    quiet = TRUE
)

# 2. Constrained Model
# Prior knowledge: residual SD should be small, around 1 g/day.
# dunif(0, 2) restricts sigma_e to [0, 2], strongly informative.
variance_prior <- list(
    sigma_e_Growth_g_day = "dunif(0, 2)"
)

fit_constrained_var <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df_small,
    priors = variance_prior,
    quiet = TRUE
)

# Visualize: The constrained posterior for sigma will be shifted left (towards 1) and sharper
plot_posterior(
    list(Default = fit_default_var, Constrained = fit_constrained_var),
    parameter = "sigma"
)

## ----check_code---------------------------------------------------------------
fit_check <- because(
    equations = list(Growth_g_day ~ Temp_Centered),
    data = df,
    n.iter = 0, quiet = TRUE
)

# Print the JAGS model
fit_check$model

