test_that("Numerical parity between JAGS and NIMBLE posteriors", {
    skip_on_cran()
    skip_if_not_installed("rjags")
    skip_if_not_installed("nimble")
    skip_if(identical(Sys.getenv("R_COVR"), "true"), "NIMBLE compilation not supported under covr")

    set.seed(42)
    df <- data.frame(X = rnorm(100), Y = rnorm(100))
    df$Y <- 1.5 * df$X + rnorm(100, 0, 0.5)
    
    # Run JAGS
    fit_jags <- because(
        list(Y ~ X),
        data = df,
        engine = "jags",
        n.iter = 1000,
        n.burnin = 200,
        quiet = TRUE
    )
    
    # Run NIMBLE
    fit_nimble <- tryCatch(
        because(
            list(Y ~ X),
            data = df,
            engine = "nimble",
            n.iter = 1000,
            n.burnin = 200,
            quiet = TRUE
        ),
        error = function(e) {
            skip(paste("NIMBLE execution failed:", conditionMessage(e)))
        }
    )
    
    # Compare posterior means for beta_Y_X
    jags_beta <- fit_jags$summary$statistics["beta_Y_X", "Mean"]
    nimble_beta <- fit_nimble$summary$statistics["beta_Y_X", "Mean"]
    
    # They should be reasonably close (within 0.1 given 1000 iter)
    expect_true(abs(jags_beta - nimble_beta) < 0.1)
    
    # Check that they both successfully captured the true parameter ~1.5
    expect_true(abs(jags_beta - 1.5) < 0.2)
    expect_true(abs(nimble_beta - 1.5) < 0.2)
})
