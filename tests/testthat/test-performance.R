test_that("Performance is acceptable for small models", {
    skip_on_cran()
    skip_if(identical(Sys.getenv("R_COVR"), "true"), "Timing benchmarks inaccurate under covr")

    set.seed(123)
    N <- 50
    X <- rnorm(N)
    Y <- 0.5 * X + rnorm(N)

    data_list <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)

    # Benchmark execution time
    time <- system.time({
        fit <- because(
            data = data_list,
            
            equations = equations,
            n.iter = 1000,
            n.burnin = 500,
            n.chains = 2,
            quiet = TRUE
        )
    })

    # Should be faster than 30 seconds for this very simple model
    # (Adjust threshold based on typical CI runner speed)
    expect_lt(time["elapsed"], 30)
})

test_that("Performance scales reasonably with N", {
    skip_on_cran()
    skip_if(identical(Sys.getenv("R_COVR"), "true"), "Timing benchmarks inaccurate under covr")

    # Compare N=20 vs N=100
    # JAGS scaling is roughly linear or N^2 depending on matrix ops

    run_model <- function(N) {
        X <- rnorm(N)
        Y <- rnorm(N)
        data <- list(X = X, Y = Y, N = N)
        equations <- list(Y ~ X)

        system.time({
            because(
                data = data,
                
                equations = equations,
                n.iter = 100,
                n.burnin = 50,
                n.chains = 2,
                quiet = TRUE
            )
        })["elapsed"]
    }

    time_20 <- run_model(20)
    time_100 <- run_model(100)

    # N=100 shouldn't be more than 30x slower than N=20 (conservative check for system variability)
    expect_lt(time_100, time_20 * 30)
})
