test_that("ic_recompile computes DIC with parallel chains", {
    skip_on_cran()

    set.seed(123)
    N <- 20
    N <- 20
    X <- rnorm(N)
    Y <- 0.5 * X + rnorm(N)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Run with parallel execution AND ic_recompile=TRUE
    fit <- because(
        data = data_list,
        
        equations = equations,
        n.iter = 500,
        n.burnin = 250,
        n.chains = 4,
        parallel = TRUE,
        n.cores = 2,
        DIC = TRUE,
        ic_recompile = TRUE,
        quiet = TRUE
    )

    # Check that DIC was computed
    expect_s3_class(fit, "because")
    expect_true(!is.null(fit$DIC))
    expect_s3_class(fit$DIC, "dic")
})

test_that("ic_recompile=FALSE skips DIC with parallel chains", {
    skip_on_cran()

    set.seed(456)
    N <- 15
    N <- 15
    X <- rnorm(N)
    Y <- rnorm(N)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Run with parallel execution but ic_recompile=FALSE
    suppressWarnings({
        fit <- because(
            data = data_list,
            
            equations = equations,
            n.iter = 200,
            n.burnin = 100,
            n.chains = 2,
            parallel = TRUE,
            n.cores = 2,
            DIC = TRUE,
            ic_recompile = FALSE,
            quiet = TRUE
        )
    })

    # DIC should be NULL
    expect_true(is.null(fit$DIC))
})
