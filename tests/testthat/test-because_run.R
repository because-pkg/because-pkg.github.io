test_that("because validates inputs", {
    # Missing data
    expect_error(
        because(data = NULL, equations = list(Y ~ X)),
        "data"
    )
})

test_that("because completes a simple run", {
    skip_on_cran() # Skip on CRAN to save time

    # Simulate simple data
    set.seed(123)
    N <- 20
    X <- rnorm(N)
    Y <- 0.5 * X + rnorm(N)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Run with very few iterations just to check it runs
    fit <- because(
        data = data_list,
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
    expect_true(!is.null(fit$model))
    expect_true(!is.null(fit$samples))
})
