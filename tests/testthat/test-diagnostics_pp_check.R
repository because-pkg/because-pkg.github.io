test_that("pp_check and bayesian diagnostics return ggplot or diagnostic objects", {
    skip_if_not_installed("bayesplot")
    skip_if_not_installed("ggplot2")

    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            c(rnorm(200, mean = 1.0, sd = 0.05), # alpha_Y
              rnorm(200, mean = 0.5, sd = 0.05), # beta_Y_X
              rep(0.5, 200)),                    # sigma_Y_res
            ncol = 3,
            dimnames = list(NULL, c("alpha_Y", "beta_Y_X", "sigma_Y_res"))
        ))
    )

    df <- data.frame(
        X = rnorm(30),
        Y = rnorm(30, mean = 1.0, sd = 0.5)
    )

    parameter_map <- data.frame(
        response = c("Y", "Y"),
        predictor = c("(Intercept)", "X"),
        parameter = c("alpha_Y", "beta_Y_X"),
        type = c("intercept", "coefficient"),
        equation_index = c(1, 1),
        stringsAsFactors = FALSE
    )

    fit <- list(
        samples = samples,
        dsep = FALSE,
        family = list(Y = "gaussian"),
        parameter_map = parameter_map,
        data = df,
        original_data = df,
        equations = list(Y ~ X)
    )
    class(fit) <- "because"

    # 1. pp_check dens_overlay
    p_dens <- pp_check(fit, resp = "Y", type = "dens_overlay", ndraws = 20)
    expect_s3_class(p_dens, "ggplot")

    # 2. pp_check hist
    p_hist <- pp_check(fit, resp = "Y", type = "hist", ndraws = 5)
    expect_s3_class(p_hist, "ggplot")

    # 3. check_normality / plot_qq
    norm_res <- check_normality(fit, resp = "Y", ndraws = 50)
    expect_s3_class(norm_res, "because_normality")
    expect_s3_class(norm_res$plot, "ggplot")
    p_qq <- plot_qq(fit, resp = "Y", ndraws = 50)
    expect_s3_class(p_qq, "ggplot")

    # 4. check_homoscedasticity
    homo_res <- check_homoscedasticity(fit, resp = "Y", ndraws = 50)
    expect_s3_class(homo_res, "because_homoscedasticity")
    expect_s3_class(homo_res$plot_residuals, "ggplot")
    expect_s3_class(homo_res$plot_binned, "ggplot")

    # 5. check_collinearity
    vif_res <- check_collinearity(fit)
    expect_true(is.data.frame(vif_res) || is.list(vif_res))
})
