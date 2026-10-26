test_that("Parallel chains run successfully", {
    skip_on_cran()
    skip_if(identical(Sys.getenv("R_COVR"), "true"), "Parallel clusters incompatible with covr")

    set.seed(123)
    N <- 20
    N <- 20
    X <- rnorm(N)
    Y <- 0.5 * X + rnorm(N)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Run with parallel execution
    fit_parallel <- because(
        data = data_list,
        
        equations = equations,
        n.iter = 200,
        n.burnin = 100,
        n.chains = 4,
        parallel = TRUE,
        n.cores = 2,
        DIC = FALSE,
        quiet = TRUE
    )

    # Check structure
    expect_s3_class(fit_parallel, "because")
    expect_s3_class(fit_parallel$samples, "mcmc.list")
    expect_equal(length(fit_parallel$samples), 4) # 4 chains
})

test_that("Sequential and parallel produce compatible results", {
    skip_on_cran()
    skip_if(identical(Sys.getenv("R_COVR"), "true"), "Parallel clusters incompatible with covr")

    set.seed(456)
    N <- 15
    N <- 15
    X <- rnorm(N)
    Y <- 0.6 * X + rnorm(N)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Sequential
    fit_seq <- because(
        data = data_list,
        
        equations = equations,
        n.iter = 200,
        n.burnin = 100,
        n.chains = 2,
        parallel = FALSE,
        quiet = TRUE
    )

    # Parallel
    fit_par <- because(
        data = data_list,
        
        equations = equations,
        n.iter = 200,
        n.burnin = 100,
        n.chains = 2,
        parallel = TRUE,
        n.cores = 2,
        DIC = FALSE,
        quiet = TRUE
    )

    # Both should produce valid results
    expect_s3_class(fit_seq, "because")
    expect_s3_class(fit_par, "because")

    # Check that parameter estimates are in similar range
    # Check that parameter estimates are in similar range
    beta_seq <- fit_seq$summary$statistics["beta_Y_X", "Mean"]
    beta_par <- fit_par$summary$statistics["beta_Y_X", "Mean"]

    # They won't be identical (different RNG streams), but should be reasonable
    expect_true(beta_seq > -0.5 && beta_seq < 2.0)
    expect_true(beta_par > -0.5 && beta_par < 2.0)
})

test_that("Parallel=FALSE with n.cores>1 still runs sequentially", {
    skip_on_cran()

    set.seed(789)
    N <- 10
    N <- 10
    X <- rnorm(N)
    Y <- rnorm(N)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # This should run sequentially even though n.cores=4
    fit <- because(
        data = data_list,
        
        equations = equations,
        n.iter = 100,
        n.burnin = 50,
        n.chains = 2,
        parallel = FALSE,
        n.cores = 4,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
})
