test_that("because_model generates valid JAGS code", {
    equations <- list(Y ~ X)

    model_output <- because_model(equations)

    # Test structure
    expect_type(model_output, "list")
    expect_named(model_output, c("model", "parameter_map"))

    # Test model string
    expect_type(model_output$model, "character")
    expect_match(model_output$model, "model \\{")
    expect_match(model_output$model, "beta_Y_X")

    # Test parameter map
    expect_s3_class(model_output$parameter_map, "data.frame")
    expect_true(nrow(model_output$parameter_map) > 0)
})

test_that("because_model handles missing data", {
    equations <- list(Y ~ X)
    vars_with_na <- "Y"

    model_output <- because_model(
        equations,
        vars_with_na = vars_with_na,
        structures = list(tree = matrix(0, 1, 1))
    )

    # Should use GLMM formulation
    expect_match(model_output$model, "tau_e_Y")
    expect_match(model_output$model, "tau_u_Y_tree")
})

test_that("because_model handles binomial variables", {
    equations <- list(Binary ~ X)
    distribution <- c(Binary = "binomial")

    model_output <- because_model(
        equations,
        family = distribution
    )

    expect_match(model_output$model, "dbern")
    expect_match(model_output$model, "logit")
})

test_that("because_model handles latent variables", {
    equations <- list(X ~ L, Y ~ L)
    induced_corr <- list(c("X", "Y"))

    model_output <- because_model(
        equations,
        induced_correlations = induced_corr
    )

    expect_match(model_output$model, "rho_X_Y")
    expect_match(model_output$model, "dmnorm")
})
