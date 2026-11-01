testthat::test_that("Negative Binomial size parameter defaults to weakly informative dgamma(2, 0.1)", {
  eqs <- list(Y ~ X)
  m <- because_model(eqs, family = list(Y = "negbinomial"))
  testthat::expect_true(any(grepl("r_Y ~ dgamma\\(2, 0.1\\)", m$model)))
})

testthat::test_that("Negative Binomial size parameter supports custom prior by parameter name", {
  eqs <- list(Y ~ X)
  m <- because_model(eqs, family = list(Y = "negbinomial"), priors = list(r_Y = "dgamma(5, 0.5)"))
  testthat::expect_true(any(grepl("r_Y ~ dgamma\\(5, 0.5\\)", m$model)))
})

testthat::test_that("Negative Binomial size parameter supports custom prior by type 'r'", {
  eqs <- list(Y ~ X)
  m <- because_model(eqs, family = list(Y = "negbinomial"), priors = list(r = "dgamma(10, 1)"))
  testthat::expect_true(any(grepl("r_Y ~ dgamma\\(10, 1\\)", m$model)))
})

testthat::test_that("Zero-inflation probability supports custom prior by parameter name and type", {
  eqs <- list(Y ~ X)
  m_def <- because_model(eqs, family = list(Y = "zinb"))
  testthat::expect_true(any(grepl("psi_Y ~ dunif\\(0, 1\\)", m_def$model)))

  m_custom <- because_model(eqs, family = list(Y = "zinb"), priors = list(psi_Y = "dbeta(2, 5)"))
  testthat::expect_true(any(grepl("psi_Y ~ dbeta\\(2, 5\\)", m_custom$model)))
})
