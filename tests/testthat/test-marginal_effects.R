test_that("marginal_effects computes AME on Gaussian and non-Gaussian models", {
    # Mock samples for a 2-equation SEM:
    # M ~ X (gaussian)
    # Y ~ M + X (poisson)
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            c(rnorm(400, mean = 0.5, sd = 0.05),
              rnorm(400, mean = 0.3, sd = 0.05),
              rnorm(400, mean = 0.2, sd = 0.05),
              rnorm(400, mean = 0.1, sd = 0.02),
              rnorm(400, mean = 0.4, sd = 0.02)),
            ncol = 5,
            dimnames = list(NULL, c("alpha_M", "beta_M_X", "alpha_Y", "beta_Y_M", "beta_Y_X"))
        ))
    )

    df <- data.frame(
        X = rnorm(30),
        M = rnorm(30),
        Y = rpois(30, lambda = 3)
    )

    parameter_map <- data.frame(
        response = c("M", "M", "Y", "Y", "Y"),
        predictor = c("(Intercept)", "X", "(Intercept)", "M", "X"),
        parameter = c("alpha_M", "beta_M_X", "alpha_Y", "beta_Y_M", "beta_Y_X"),
        type = c("intercept", "coefficient", "intercept", "coefficient", "coefficient"),
        equation_index = c(1, 1, 2, 2, 2),
        stringsAsFactors = FALSE
    )

    fit <- list(
        samples = samples,
        dsep = FALSE,
        family = list(M = "gaussian", Y = "poisson"),
        parameter_map = parameter_map,
        data = df,
        original_data = df,
        equations = list(M ~ X, Y ~ M + X)
    )
    class(fit) <- "because"

    # 1. Compute Average Marginal Effects (default)
    me <- marginal_effects(fit, prob = 0.95)

    expect_s3_class(me, "data.frame")
    expect_named(me, c("Response", "Predictor", "Category", "Effect", "Lower", "Upper", "Family"))
    expect_equal(nrow(me), 3) # beta_M_X, beta_Y_M, beta_Y_X
    expect_equal(me$Response, c("M", "Y", "Y"))
    expect_equal(me$Predictor, c("X", "M", "X"))
    expect_true(all(is.na(me$Category))) # Continuous variables have NA category
    expect_equal(me$Family, c("gaussian", "poisson", "poisson"))

    # For Gaussian M ~ X: AME equals the linear slope ~0.3
    expect_equal(me$Effect[me$Response == "M"], 0.3, tolerance = 0.05)

    # For Poisson Y: AME is on expected count scale (positive and proportional to count)
    expect_true(me$Effect[me$Response == "Y" & me$Predictor == "M"] > 0)

    # 2. Compute Marginal Effects at the Mean (MEM)
    me_mem <- marginal_effects(fit, at = "mean")
    expect_s3_class(me_mem, "data.frame")
    expect_equal(nrow(me_mem), 3)

    # 3. Custom credible interval probability
    me_90 <- marginal_effects(fit, prob = 0.90)
    expect_true(all(me_90$Lower >= me$Lower))
    expect_true(all(me_90$Upper <= me$Upper))
})

test_that("marginal_effects handles categorical factor predictors", {
    # Predictor with 2 categories (Treatment vs Control)
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            c(rnorm(200, mean = 0.5), rnorm(200, mean = 1.2)),
            ncol = 2,
            dimnames = list(NULL, c("alpha_Y", "beta_Y_Trt_Treated"))
        ))
    )

    df <- data.frame(
        Trt = factor(rep(c("Control", "Treated"), each = 15)),
        Trt_Treated = rep(c(0, 1), each = 15),
        Y = rnorm(30)
    )

    parameter_map <- data.frame(
        response = c("Y", "Y"),
        predictor = c("(Intercept)", "Trt_Treated"),
        parameter = c("alpha_Y", "beta_Y_Trt_Treated"),
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
        equations = list(Y ~ Trt),
        categorical_vars = list(Trt = list(type = "unordered", levels = c("Control", "Treated")))
    )
    class(fit) <- "because"

    me <- marginal_effects(fit)
    expect_s3_class(me, "data.frame")
    expect_equal(me$Predictor, "Trt")
    expect_equal(me$Category, "Treated")
    expect_equal(me$Effect, 1.2, tolerance = 0.05)
})
