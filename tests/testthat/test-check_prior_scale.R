testthat::test_that("check_prior_scale warns when continuous traits have large SD", {
  df <- data.frame(
    Y = rnorm(50, mean = 100, sd = 50),
    X = rnorm(50, mean = 0, sd = 1)
  )
  eqs <- list(Y ~ X)

  testthat::expect_warning(
    check_prior_scale(df, eqs, family = list(Y = "gaussian"), priors = NULL, quiet = FALSE),
    "Continuous variable\\(s\\) with large variance detected: 'Y'"
  )
})

testthat::test_that("check_prior_scale does not warn for standardized continuous traits", {
  df <- data.frame(
    Y = rnorm(50, mean = 0, sd = 1),
    X = rnorm(50, mean = 0, sd = 1)
  )
  eqs <- list(Y ~ X)

  testthat::expect_silent(
    check_prior_scale(df, eqs, family = list(Y = "gaussian"), priors = NULL, quiet = FALSE)
  )
})

testthat::test_that("check_prior_scale does not warn when custom priors are provided", {
  df <- data.frame(
    Y = rnorm(50, mean = 100, sd = 50),
    X = rnorm(50, mean = 0, sd = 1)
  )
  eqs <- list(Y ~ X)
  custom_priors <- list(sigma_Y_res = "dunif(0, 500)")

  testthat::expect_silent(
    check_prior_scale(df, eqs, family = list(Y = "gaussian"), priors = custom_priors, quiet = FALSE)
  )
})

testthat::test_that("check_prior_scale does not warn for non-Gaussian (discrete/count) traits", {
  df <- data.frame(
    Y = rpois(50, lambda = 100),
    X = rnorm(50, mean = 0, sd = 1)
  )
  # Even if count SD > 20
  eqs <- list(Y ~ X)

  testthat::expect_silent(
    check_prior_scale(df, eqs, family = list(Y = "poisson"), priors = NULL, quiet = FALSE)
  )
})

testthat::test_that("check_prior_scale respects quiet = TRUE", {
  df <- data.frame(
    Y = rnorm(50, mean = 100, sd = 50),
    X = rnorm(50, mean = 0, sd = 1)
  )
  eqs <- list(Y ~ X)

  testthat::expect_silent(
    check_prior_scale(df, eqs, family = list(Y = "gaussian"), priors = NULL, quiet = TRUE)
  )
})
