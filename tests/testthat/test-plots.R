test_that("plot_coef renders ggplot objects for raw and marginal", {
    skip_if_not_installed("ggplot2")

    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            c(rnorm(100, mean = 0.5, sd = 0.05), # beta_Y_X
              rnorm(100, mean = 0.0, sd = 0.05)),# beta_Y_Z
            ncol = 2,
            dimnames = list(NULL, c("beta_Y_X", "beta_Y_Z"))
        ))
    )

    df <- data.frame(
        X = rnorm(30),
        Z = rnorm(30),
        Y = rnorm(30)
    )

    parameter_map <- data.frame(
        response = c("Y", "Y"),
        predictor = c("X", "Z"),
        parameter = c("beta_Y_X", "beta_Y_Z"),
        type = c("coefficient", "coefficient"),
        equation_index = c(1, 1),
        stringsAsFactors = FALSE
    )

    summ <- summary(samples)

    fit <- list(
        samples = samples,
        summary = summ,
        dsep = FALSE,
        family = list(Y = "gaussian"),
        parameter_map = parameter_map,
        data = df,
        original_data = df,
        equations = list(Y ~ X + Z)
    )
    class(fit) <- "because"

    # Raw plot
    p_raw <- plot_coef(fit, type = "raw")
    expect_s3_class(p_raw, "ggplot")

    # Marginal plot
    p_mar <- plot_coef(fit, type = "marginal")
    expect_s3_class(p_mar, "ggplot")
})

test_that("plot_posterior executes and generates plots", {
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(100, mean = 0.5, sd = 0.05),
            ncol = 1,
            dimnames = list(NULL, c("beta_Y_X"))
        ))
    )

    fit <- list(
        samples = samples,
        equations = list(Y ~ X)
    )
    class(fit) <- "because"

    # Should run silently without error
    expect_error(plot_posterior(fit, parameter = "beta_Y_X"), NA)
})

test_that("plot_dag generates a ggplot visualization", {
    skip_if_not_installed("dagitty")
    skip_if_not_installed("ggdag")
    skip_if_not_installed("ggraph")
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("dplyr")

    eqs <- list(Y ~ X + Z, X ~ Z)
    p_dag <- plot_dag(eqs)
    expect_s3_class(p_dag, "ggplot")
})
