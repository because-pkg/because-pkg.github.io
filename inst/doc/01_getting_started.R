## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", eval = FALSE
)

## -----------------------------------------------------------------------------
# remotes::install_github("because-pkg/because@v1.2.4", build_vignettes = TRUE)

## -----------------------------------------------------------------------------
# remotes::install_github("because-pkg/because", build_vignettes = TRUE)

## -----------------------------------------------------------------------------
# library(because)

## -----------------------------------------------------------------------------
# # set seed for reproducibility
# set.seed(67)
# 
# # Simulate predictor
# X <- rnorm(n = 100, mean = 50, sd = 10)
# 
# # Generate response with the chosen intercept (alpha) and slope (beta)
# alpha <- 20
# Y <- alpha + beta * X + rnorm(n, mean = 0, sd = 10)
# 
# # Combine into data frame
# sim.dat <- data.frame(X, Y)
# 
# # Fit linear model with lm() function for comparison
# summary(lm(Y ~ X, data = sim.dat))

## -----------------------------------------------------------------------------
# # Define the equations
# equations <- list(Y ~ X)
# 
# # Fit the model
# fit <- because(
#   equations = equations,
#   data = sim.dat
# )
# 
# # View results
# summary(fit)

## -----------------------------------------------------------------------------
# # Plot trace plots
# plot(fit$samples)

## ----echo=FALSE, eval=TRUE, out.width="100%"----------------------------------
knitr::include_graphics("figures/traceplot_1.png")

## -----------------------------------------------------------------------------
# # Generate JAGS model code
# jags_model_code <- because_model(
#   equations = equations
# )
# 
# cat(jags_model_code$model)

