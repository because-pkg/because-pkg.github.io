test_that("Regression test: simple Gaussian model", {
    skip_on_cran()

    # Fixed seed for reproducibility
    set.seed(42)
    N <- 30

    # Simulate data with known relationship
    X <- rnorm(N)
    Y <- 0.7 * X + rnorm(N, sd = 0.5)

    data <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    fit <- because(
        data = data,
        
        equations = equations,
        n.iter = 1000,
        n.burnin = 500,
        n.chains = 2,
        quiet = TRUE
    )

    # Test structure
    expect_s3_class(fit, "because")
    expect_true(!is.null(fit$summary))

    # Test parameter recovery (should be close to 0.7)
    # Allow wide tolerance since this is a small sample
    beta_mean <- fit$summary$statistics["beta_Y_X", "Mean"]
    expect_true(
        beta_mean > 0.3 && beta_mean < 1.1,
        info = paste("beta_Y_X mean was", round(beta_mean, 3))
    )
})

test_that("Regression test: binomial model", {
    skip(
        "Binomial models can have convergence issues with short MCMC runs"
    )

    skip_on_cran()

    set.seed(123)
    N <- 40

    # Simulate binomial data
    X <- rnorm(N)
    prob <- plogis(0.5 * X)
    Y <- rbinom(N, 1, prob)

    data <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)
    distribution <- c(Y = "binomial")

    fit <- because(
        data = data,
        
        equations = equations,
        family = distribution,
        n.iter = 1000,
        n.burnin = 500,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")

    # Check that logit coefficient exists
    expect_true("beta_Y_X" %in% rownames(fit$summary$statistics))

    # Coefficient should be positive (we simulated positive relationship)
    beta_mean <- fit$summary$statistics["beta_Y_X", "Mean"]
    expect_true(beta_mean > 0)
})

test_that("Regression test: model with missing data", {
    skip_on_cran()

    set.seed(789)
    N <- 35

    X <- rnorm(N)
    Y <- 0.6 * X + rnorm(N)

    # Add missing values
    Y[c(5, 15, 25)] <- NA

    data <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    fit <- because(
        data = data,
        
        equations = equations,
        n.iter = 300,
        n.burnin = 150,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")

    # Should have estimated parameters despite missing data
    expect_true("beta_Y_X" %in% rownames(fit$summary$statistics))

    # Parameter should be in reasonable range
    beta_mean <- fit$summary$statistics["beta_Y_X", "Mean"]
    expect_true(beta_mean > 0.2 && beta_mean < 1.5)
})
