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

testthat::test_that("Multinomial and ordinal priors default to regularized link scales and honor custom overrides", {
  eqs <- list(Y ~ X)
  dummy_struct <- diag(5)

  m_multi <- because_model(
    eqs,
    family = list(Y = "multinomial"),
    structures = list(spatial = dummy_struct)
  )
  testthat::expect_true(any(grepl("alpha_Y\\[k\\] ~ dnorm\\(0, 1\\)", m_multi$model)))
  testthat::expect_true(any(grepl("beta_Y_X\\[k\\] ~ dnorm\\(0, 1\\)", m_multi$model)))
  testthat::expect_true(any(grepl("tau_u_spatial_Y\\[k\\] ~ dgamma\\(10, 10\\)", m_multi$model)))
  testthat::expect_false(any(grepl("dgamma\\(0\\.001, 0\\.001\\)", m_multi$model)))

  m_multi_custom <- because_model(
    eqs,
    family = list(Y = "multinomial"),
    structures = list(spatial = dummy_struct),
    priors = list(
      alpha_Y = "dnorm(0, 0.5)",
      beta_Y_X = "dnorm(0, 0.25)",
      tau_u_spatial_Y = "dgamma(2, 2)"
    )
  )
  testthat::expect_true(any(grepl("alpha_Y\\[k\\] ~ dnorm\\(0, 0\\.5\\)", m_multi_custom$model)))
  testthat::expect_true(any(grepl("beta_Y_X\\[k\\] ~ dnorm\\(0, 0\\.25\\)", m_multi_custom$model)))
  testthat::expect_true(any(grepl("tau_u_spatial_Y\\[k\\] ~ dgamma\\(2, 2\\)", m_multi_custom$model)))

  m_ord <- because_model(
    eqs,
    family = list(Y = "ordinal"),
    structures = list(spatial = dummy_struct)
  )
  testthat::expect_true(any(grepl("beta_Y_X ~ dnorm\\(0, 1\\)", m_ord$model)))
  testthat::expect_true(any(grepl("tau_u_spatial_Y ~ dgamma\\(10, 10\\)", m_ord$model)))
  testthat::expect_false(any(grepl("dgamma\\(0\\.001, 0\\.001\\)", m_ord$model)))

  m_ord_custom <- because_model(
    eqs,
    family = list(Y = "ordinal"),
    structures = list(spatial = dummy_struct),
    priors = list(
      cutpoint = "dnorm(0, 0.5)",
      beta_Y_X = "dnorm(0, 0.25)",
      tau_u_spatial_Y = "dgamma(2, 2)"
    )
  )
  testthat::expect_true(any(grepl("cutpoint_raw_Y\\[1\\] ~ dnorm\\(0, 0\\.5\\)", m_ord_custom$model)))
  testthat::expect_true(any(grepl("beta_Y_X ~ dnorm\\(0, 0\\.25\\)", m_ord_custom$model)))
  testthat::expect_true(any(grepl("tau_u_spatial_Y ~ dgamma\\(2, 2\\)", m_ord_custom$model)))
})
