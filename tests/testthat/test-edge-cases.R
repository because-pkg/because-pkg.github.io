test_that("Handles small sample size gracefully", {
    # N=2 which is the minimum for meaningful statistics
    data <- list(X = c(1,2), Y = c(1,2), N = 2)
    equations <- list(Y ~ X)

    fit <- because(
        data = data,
        equations = equations,
        n.chains = 2,
        n.iter = 100, n.thin = 1,
        n.burnin = 0,
        quiet = TRUE
    )
    expect_s3_class(fit, "because")
})

test_that("Handles perfect collinearity", {
    skip_on_cran()

    set.seed(123)
    N <- 20
    N <- 20
    X <- rnorm(N)
    Y <- 2 * X # Perfect correlation

    data <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Should run but might have convergence issues or high variance
    # We just want to ensure it doesn't crash
    fit <- because(
        data = data,
        
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
})

test_that("Handles disconnected DAG", {
    # X -> Y, Z (Z is independent)
    set.seed(123)
    N <- 20
    X <- rnorm(N)
    Y <- X + rnorm(N)
    Z <- rnorm(N)

    data <- list(X = X, Y = Y, Z = Z, N = N)
    equations <- list(
        Y ~ X,
        Z ~ 1 # Intercept only
    )

    fit <- because(
        data = data,
        
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
})
