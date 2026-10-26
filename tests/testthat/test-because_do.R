test_that("because_do simulates downstream interventions correctly", {
    set.seed(123)
    # 2-step causal mediation chain:
    # X -> M -> Y
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            c(rnorm(400, mean = 0.0), # alpha_M
              rnorm(400, mean = 0.5), # beta_M_X
              rnorm(400, mean = 0.0), # alpha_Y
              rnorm(400, mean = 0.8), # beta_Y_M
              rep(0.01, 400),         # sigma_M
              rep(0.01, 400)),        # sigma_Y
            ncol = 6,
            dimnames = list(NULL, c("alpha_M", "beta_M_X", "alpha_Y", "beta_Y_M", "sigma_M_res", "sigma_Y_res"))
        ))
    )

    df <- data.frame(
        X = rnorm(50, mean = 0, sd = 1),
        M = rnorm(50),
        Y = rnorm(50)
    )

    parameter_map <- data.frame(
        response = c("M", "M", "Y", "Y"),
        predictor = c("(Intercept)", "X", "(Intercept)", "M"),
        parameter = c("alpha_M", "beta_M_X", "alpha_Y", "beta_Y_M"),
        type = c("intercept", "coefficient", "intercept", "coefficient"),
        equation_index = c(1, 1, 2, 2),
        stringsAsFactors = FALSE
    )

    fit <- list(
        samples = samples,
        dsep = FALSE,
        family = list(M = "gaussian", Y = "gaussian"),
        parameter_map = parameter_map,
        data = df,
        original_data = df,
        equations = list(M ~ X, Y ~ M)
    )
    class(fit) <- "because"

    # 1. Formula intervention: shift X by +1 SD
    cf_formula <- because_do(fit, X = ~ . + 1)
    expect_s3_class(cf_formula, "because_counterfactual")
    expect_true(is.matrix(cf_formula$X))

    summ_formula <- summary(cf_formula)
    expect_s3_class(summ_formula, "data.frame")
    expect_true(all(c("Variable", "Mean", "SD", "2.5%", "50%", "97.5%") %in% colnames(summ_formula)))

    # When X is shifted by +1:
    # M shifts by ~ +0.5
    # Y shifts by ~ 0.5 * 0.8 = +0.4
    expect_equal(summ_formula$Mean[summ_formula$Variable == "X"], mean(df$X) + 1, tolerance = 0.05)
    expect_equal(summ_formula$Mean[summ_formula$Variable == "M"], 0.5 * (mean(df$X) + 1), tolerance = 0.1)

    # 2. Direct atomic setting (clamp X = 2)
    cf_atomic <- because_do(fit, X = 2)
    summ_atomic <- summary(cf_atomic)
    expect_equal(summ_atomic$Mean[summ_atomic$Variable == "X"], 2.0, tolerance = 0.001)
    expect_equal(summ_atomic$SD[summ_atomic$Variable == "X"], 0.0, tolerance = 0.001)

    # 3. List wrapper syntax
    cf_list <- because_do(fit, interventions = list(X = quote(X + 1)))
    expect_s3_class(cf_list, "because_counterfactual")
    summ_list <- summary(cf_list)
    expect_equal(summ_list$Mean[summ_list$Variable == "X"], mean(df$X) + 1, tolerance = 0.05)
})
