library(testthat)
library(because)

test_that("because() calculates Rhat and handles non-convergence", {
    set.seed(42)
    df <- data.frame(X = rnorm(50), Y = rnorm(50))
    
    # Run with very few iterations to force poor convergence
    # We use quiet=FALSE to ensure the warnings are printed
    expect_warning(
        fit <- because(
            list(Y ~ X),
            data = df,
            engine = "jags",
            n.iter = 15,
            n.burnin = 2,
            n.chains = 2,
            quiet = FALSE
        ),
        regexp = "Adaptation incomplete|R-hat"
    )
    
    expect_s3_class(fit, "because")
})
