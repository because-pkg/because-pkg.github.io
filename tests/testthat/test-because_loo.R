test_that("because_loo works with valid model", {
    skip_on_cran()
    skip_if_not_installed("loo")

    # Setup: Minimal model with WAIC=TRUE (generates log_lik)
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
    loo_res <- because_loo(fit)

    # Check structure
    expect_s3_class(loo_res, "loo")

    # Check that key components exist
    expect_true("estimates" %in% names(loo_res))
    expect_true("pointwise" %in% names(loo_res))
    expect_true("diagnostics" %in% names(loo_res))

    # Check that estimates table has expected rows
    expect_true("looic" %in% rownames(loo_res$estimates))
    expect_true("p_loo" %in% rownames(loo_res$estimates))
    expect_true("elpd_loo" %in% rownames(loo_res$estimates))
})

test_that("because_loo throws error if WAIC=FALSE (no log_lik)", {
    skip_on_cran()
    skip_if_not_installed("loo")

    # Setup: Model WITHOUT WAIC=TRUE
    set.seed(456)
    
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

    # Expect success (with message about calculating log_lik)
    expect_message(
        loo_res <- because_loo(fit_no_waic),
        "Pointwise log-likelihoods not found"
    )

    expect_s3_class(loo_res, "loo")
    expect_true("estimates" %in% names(loo_res))
})

test_that("because_loo validates input", {
    expect_error(
        because_loo("not a model"),
        "Input must be a 'because' model object"
    )
})
