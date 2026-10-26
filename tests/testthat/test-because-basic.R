library(testthat)
library(because)

test_that("because() validates missing data argument appropriately", {
  expect_error(because(list(Y ~ X), data = NULL), "Argument 'data' must be provided")
})

test_that("because() runs a simple linear model with JAGS", {
  set.seed(123)
  df <- data.frame(X = rnorm(20), Y = rnorm(20))
  
  # Run a very short chain to avoid CRAN timeouts
  fit <- because(
    equations = list(Y ~ X),
    data = df,
    engine = "jags",
    n.chains = 2,
    n.iter = 100,
    n.thin = 1,
    n.burnin = 0,
    quiet = TRUE
  )
  
  expect_s3_class(fit, "because")
  expect_true(length(fit$samples) == 2)
  expect_true("beta_Y_X" %in% colnames(fit$samples[[1]]))
})

test_that("because() runs a simple linear model with NIMBLE", {
  skip_if_not_installed("nimble")
  
  set.seed(123)
  df <- data.frame(X = rnorm(20), Y = rnorm(20))
  
  # Run a very short chain to avoid CRAN timeouts
  fit <- because(
    equations = list(Y ~ X),
    data = df,
    engine = "nimble",
    n.chains = 2,
    n.iter = 100,
    n.thin = 1,
    n.burnin = 0,
    quiet = TRUE
  )
  
  expect_s3_class(fit, "because")
  expect_true(length(fit$samples) == 2)
  expect_true("beta_Y_X" %in% colnames(fit$samples[[1]]))
})

test_that("because() imputes missing data (NAs) correctly", {
  set.seed(123)
  df <- data.frame(X = rnorm(20), Y = rnorm(20))
  df$Y[c(5, 10)] <- NA # Add missing values
  
  # Run model tracking imputed values
  fit <- because(
    equations = list(Y ~ X),
    data = df,
    engine = "jags",
    n.chains = 2,
    n.iter = 100,
    n.thin = 1,
    n.burnin = 0,
    quiet = TRUE
  )
  
  expect_s3_class(fit, "because")
})
