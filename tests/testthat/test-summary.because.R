test_that("summary.because handles standard output", {
    # Mock a because object
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(100),
            ncol = 2,
            dimnames = list(NULL, c("alpha", "beta"))
        )),
        coda::mcmc(matrix(
            rnorm(100),
            ncol = 2,
            dimnames = list(NULL, c("alpha", "beta"))
        ))
    )

    fit <- list(
        samples = samples,
        dsep = FALSE,
        DIC = c(deviance = 100, pD = 2, DIC = 102),
        WAIC = c(waic = 100, p_waic = 2)
    )
    class(fit) <- "because"

    # Capture output to verify printing (explicit print required)
    summ <- summary(fit)
    output <- capture_output({
        print(summ)
    })

    # summary.because returns a summary.because object
    expect_s3_class(summ, "summary.because")
    expect_true(is.matrix(summ$results) || is.data.frame(summ$results))
    expect_match(output, "DIC")
    expect_match(output, "WAIC")
})

test_that("summary.because handles d-sep output", {
    # Mock a d-sep because object
    # Use 2 parameters to ensure summary returns a matrix
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(400, mean = 0),
            ncol = 2,
            dimnames = list(NULL, c("beta_Y_X", "beta_Z_Y"))
        ))
    )

    # Setup d-sep specific fields
    dsep_tests <- list(Y ~ X)
    attr(dsep_tests[[1]], "test_var") <- "X"

    parameter_map <- data.frame(
        response = "Y",
        predictor = "X",
        parameter = "beta_Y_X",
        equation_index = 1,
        stringsAsFactors = FALSE
    )

    fit <- list(
        samples = samples,
        dsep = TRUE,
        dsep_tests = dsep_tests,
        parameter_map = parameter_map,
        dsep_results = list(
            list(
                samples = samples,
                param_map = parameter_map
            )
        )
    )
    class(fit) <- "because"

    results <- summary(fit)
    output <- capture_output({
        print(results)
    })

    # Output matches "d-separation Tests" (updated string)
    expect_match(output, "d-separation Tests")
    expect_match(output, "Y _\\|\\|_ X")

    # Check structure of returned object
    expect_s3_class(results, "summary.because")
    expect_true(is.data.frame(results$results))
    expect_true(nrow(results$results) == 1)
    expect_equal(ncol(results$results), 7)
    expect_true("LowerCI" %in% colnames(results$results))
})
