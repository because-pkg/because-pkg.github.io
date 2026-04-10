## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
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

## ----standardize--------------------------------------------------------------
# Standardize all variables
eco_data_std <- scale(eco_data)
head(eco_data_std)

## ----fit_model----------------------------------------------------------------
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

## ----plot_dag-----------------------------------------------------------------
plot_dag(fit)

## ----mediation----------------------------------------------------------------
# Run Mediation Analysis for Elevation -> Abundance
med_results <- because_mediation(fit, exposure = "Elevation", outcome = "Abundance")

## ----summary------------------------------------------------------------------
med_results$summary

## ----paths--------------------------------------------------------------------
med_results$paths

