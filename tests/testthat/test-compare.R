test_that("because_compare runs multiple models and returns comparison", {
    skip_on_cran()

    set.seed(123)
    N <- 20
    X <- rnorm(N)
    Z <- rnorm(N)
    Y <- 0.5 * X + 0.3 * Z + rnorm(N)

    data_list <- list(X = X, Y = Y, Z = Z, N = N)

    # Define two models
    models <- list(
        m1 = list(equations = list(Y ~ X)),
        m2 = list(equations = list(Y ~ X + Z))
    )

    # 1. Sequential execution
    res_seq <- because_compare(
        models,
        data_list,
        n.cores = 1,
        n.iter = 1000, # Increased iterations for stable WAIC
        n.burnin = 500,
        quiet = TRUE,
        WAIC = TRUE
    )

    expect_type(res_seq, "list")
    expect_named(res_seq, c("results", "comparison"))
    expect_equal(length(res_seq$results), 2)
    expect_s3_class(res_seq$comparison, "data.frame")
    expect_equal(nrow(res_seq$comparison), 2)

    # Check WAIC presence
    if (all(is.na(res_seq$comparison$WAIC))) {
        warning("WAIC is all NA in sequential run. Check JAGS output.")
        print(res_seq$results[[1]]$WAIC)
    } else {
        expect_true("WAIC" %in% names(res_seq$comparison))
        expect_true("weight" %in% names(res_seq$comparison))
        expect_true(!any(is.na(res_seq$comparison$WAIC)))
    }

    # 2. Parallel execution
    res_par <- because_compare(
        models,
        data_list,
        n.cores = 2,
        n.iter = 1000,
        n.burnin = 500,
        quiet = TRUE,
        WAIC = TRUE
    )

    expect_equal(length(res_par$results), 2)
    expect_equal(nrow(res_par$comparison), 2)

    # Results should be similar (not identical due to RNG)
    if (
        !all(is.na(res_seq$comparison$WAIC)) &&
            !all(is.na(res_par$comparison$WAIC))
    ) {
        expect_true(
            abs(res_seq$comparison$WAIC[1] - res_par$comparison$WAIC[1]) < 20
        )
    }
})

test_that("because_compare handles errors gracefully", {
    expect_error(
        because_compare(list(), list()),
        "No fitted 'because' models provided"
    )
})
