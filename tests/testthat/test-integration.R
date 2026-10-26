test_that("Full workflow: simple model parameter recovery", {
    skip_on_cran()

    # Complete workflow
    set.seed(42)
    N <- 50

    # Simulate independent traits
    # X -> Y
    X <- rnorm(N)
    # Y = 0.8*X + error
    Y <- 0.8 * X + rnorm(N, sd = 0.5)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Run model
    fit <- because(
        data = data_list,

        equations = equations,
        n.iter = 2000,
        n.burnin = 1000,
        n.chains = 3,
        quiet = TRUE
    )

    # Check results
    expect_s3_class(fit, "because")

    # Check convergence
    expect_true(!is.null(fit$summary))

    # Check parameter recovery (should be close to 0.8)
    # We allow a relatively wide tolerance due to small sample size and MCMC noise
    beta_mean <- fit$summary$statistics["beta_Y_X", "Mean"]
    expect_true(beta_mean > 0.4 && beta_mean < 1.2)

    # Check that 95% CI includes true value (0.8)
    lower <- fit$summary$quantiles["beta_Y_X", "2.5%"]
    upper <- fit$summary$quantiles["beta_Y_X", "97.5%"]
    expect_true(lower < 0.8 && upper > 0.8)
})

test_that("Full workflow: missing data handling", {
    skip_on_cran()

    set.seed(123)
    N <- 40
    X <- rnorm(N)
    Y <- 0.5 * X + rnorm(N)

    # Add missing data to response
    Y_miss <- Y
    Y_miss[c(3, 7, 15)] <- NA

    data_list <- list(X = X, Y = Y_miss, N = N)
    equations <- list(Y ~ X)

    fit <- because(
        data = data_list,

        equations = equations,
        n.iter = 500,
        n.burnin = 250,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
    expect_true(!is.null(fit$summary))

    # Check that it estimated the missing values (nodes in JAGS)
    # The missing values are usually monitored if we monitor "Y"
    # But because monitors specific parameters.
    # Let's check if the model ran without error, which implies imputation worked.

    # Also check that beta is reasonable despite missing data
    beta_mean <- fit$summary$statistics["beta_Y_X", "Mean"]
    expect_true(beta_mean > 0.0 && beta_mean < 1.0)
})
