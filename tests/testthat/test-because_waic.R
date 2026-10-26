test_that("because_waic calculates WAIC correctly", {
    skip_on_cran()

    # Setup: Run a minimal model to get a valid because object with an rjags model
    set.seed(123)
    X <- rnorm(10)
    Y <- 0.5 * X + rnorm(10)
    data <- list(X = X, Y = Y, N = 10)
    equations <- list(Y ~ X)

    fit <- because(
        data = data,
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        quiet = TRUE,
        WAIC = TRUE
    )

    # Test execution
    waic_res <- because_waic(fit)

    # Check structure
    expect_s3_class(waic_res, "because_waic")
    expect_named(waic_res, c("Estimate", "SE"))
    expect_true(!is.na(waic_res["waic", "Estimate"]))
    expect_true(!is.na(waic_res["p_waic", "Estimate"]))
})

test_that("because_waic auto-refits when WAIC=FALSE", {
    skip_on_cran()

    set.seed(123)
    X <- rnorm(10)
    Y <- 0.5 * X + rnorm(10)
    data <- list(X = X, Y = Y, N = 10)
    equations <- list(Y ~ X)

    fit_no_waic <- because(
        data = data,
        equations = equations,
        n.iter = 50,
        n.burnin = 10,
        n.chains = 2,
        quiet = TRUE,
        WAIC = FALSE
    )

    expect_message(
        waic_res <- because_waic(fit_no_waic),
        "Pointwise log-likelihoods not found"
    )

    expect_s3_class(waic_res, "data.frame")
    expect_true(!is.na(waic_res["waic", "Estimate"]))
})

test_that("because_waic validates input", {
    expect_error(
        because_waic("not a model"),
        "Input must be a 'because' model object"
    )
})
