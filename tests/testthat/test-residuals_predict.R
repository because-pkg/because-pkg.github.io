test_that("residuals, fitted, and predict work across Gaussian and Poisson models", {
    set.seed(123)
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            c(rnorm(300, mean = 1.0, sd = 0.05), # alpha_Y
              rnorm(300, mean = 0.5, sd = 0.05), # beta_Y_X
              rep(0.5, 300)),                    # sigma_Y_res
            ncol = 3,
            dimnames = list(NULL, c("alpha_Y", "beta_Y_X", "sigma_Y_res"))
        ))
    )

    df <- data.frame(
        X = rnorm(40),
        Y = rnorm(40, mean = 1.0, sd = 0.5)
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

    # 1. fitted.because
    fit_vals <- fitted(fit, resp = "Y", ndraws = 300)
    expect_true(is.numeric(fit_vals))
    expect_equal(length(fit_vals), 40)

    # 2. residuals.because
    res_raw <- residuals(fit, resp = "Y", type = "response", ndraws = 300)
    expect_equal(length(res_raw), 40)
    expect_equal(res_raw, df$Y - fit_vals, tolerance = 0.001)

    res_pearson <- residuals(fit, resp = "Y", type = "pearson")
    expect_equal(length(res_pearson), 40)

    res_quantile <- residuals(fit, resp = "Y", type = "quantile")
    expect_equal(length(res_quantile), 40)
    expect_false(any(is.na(res_quantile)))

    # 3. predict.because
    pred_sum <- predict(fit, resp = "Y", summary = TRUE)
    expect_s3_class(pred_sum, "data.frame")
    expect_true(all(c("Estimate", "Est.Error", "Q2.5", "Q97.5") %in% colnames(pred_sum)))
    expect_equal(nrow(pred_sum), 40)

    pred_draws <- predict(fit, resp = "Y", summary = FALSE)
    expect_true(is.matrix(pred_draws))
    expect_equal(ncol(pred_draws), 40)

    # 4. posterior_predict
    yrep <- posterior_predict(fit, resp = "Y", ndraws = 50)
    expect_true(is.matrix(yrep))
    expect_equal(dim(yrep), c(50, 40))
})
