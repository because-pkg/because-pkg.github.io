test_that("plot_dsep returns a ggplot object", {
    # Mock a d-sep because object similar to test-summary.because.R
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
    attr(dsep_tests[[1]], "formula_string") <- "Y _||_ X"

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

    # plot_dsep should work on this mock object
    p <- plot_dsep(fit)

    expect_s3_class(p, "ggplot")
    expect_equal(p$labels$title, "d-separation Independence Tests")

    # Check if the data in plot has our test
    expect_true(any(grepl("Y _\\|\\|_ X", p$data$Test)))
})
