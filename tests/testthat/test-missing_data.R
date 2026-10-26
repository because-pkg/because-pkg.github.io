library(testthat)
library(because)

test_that("because() handles endogenous missing data imputation", {
    set.seed(123)
    df <- data.frame(X = rnorm(30), Y = rnorm(30))
    # Make Y endogenous with missing values
    df$Y[c(5, 15, 25)] <- NA
    
    fit <- because(
        list(Y ~ X),
        data = df,
        engine = "jags",
        n.iter = 100,
        n.burnin = 10,
        quiet = TRUE
    )
    
    expect_s3_class(fit, "because")
    
    # Check if we can extract imputed values
    # NOTE: by default, imputed values are only tracked if monitor="all" or specifically requested.
    # We will test extract_imputed with a specific run that monitors 'Y'
    fit_tracked <- because(
        list(Y ~ X),
        data = df,
        engine = "jags",
        n.iter = 100,
        n.burnin = 10,
        monitor = "all",
        quiet = TRUE
    )
    
    imputed <- extract_imputed(fit_tracked)
    expect_s3_class(imputed, "data.frame")
    expect_true(nrow(imputed) == 3)
})
