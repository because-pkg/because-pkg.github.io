test_that("Gaussian distribution generates correct JAGS code", {
    equations <- list(Y ~ X)

    model_output <- because_model(equations)

    # Should use dnorm (normal distribution)
    expect_match(model_output$model, "dnorm")

    # Should NOT use logit or dbern
    expect_false(grepl("logit", model_output$model))
    expect_false(grepl("dbern", model_output$model))

    # Should have standard beta parameters
    expect_match(model_output$model, "beta_Y_X")
})

test_that("Binomial distribution generates correct JAGS code", {
    equations <- list(BinaryOutcome ~ X)
    distribution <- c(BinaryOutcome = "binomial")

    model_output <- because_model(
        equations,
        family = distribution
    )

    # Should use dbern (Bernoulli distribution)
    expect_match(model_output$model, "dbern")

    # Should use logit link function
    expect_match(model_output$model, "logit")

    # Should have beta parameters
    expect_match(model_output$model, "beta_BinaryOutcome_X")
})

test_that("Mixed distributions work correctly", {
    # Test that we can have both Gaussian and binomial in same model
    equations <- list(
        Binary ~ X
    )
    distribution <- c(Binary = "binomial")

    model_output <- because_model(
        equations,
        family = distribution
    )

    # Should have dbern for binomial
    expect_match(model_output$model, "dbern")
    expect_match(model_output$model, "logit")
})

test_that("Binomial distribution runs successfully", {
    skip_on_cran()

    set.seed(123)
    N <- 30
    X <- rnorm(N)

    # Simulate binary outcome
    prob <- plogis(0.5 * X)
    Y <- rbinom(N, 1, prob)

    data <- list(X = X, Y = Y, N = N)
    equations <- list(Y ~ X)
    distribution <- c(Y = "binomial")

    fit <- because(
        data = data,
        
        equations = equations,
        family = distribution,
        n.iter = 200,
        n.burnin = 100,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
    expect_true("beta_Y_X" %in% rownames(fit$summary$statistics))
})

test_that("Default distribution is Gaussian", {
    equations <- list(Y ~ X)

    # Not specifying distribution should default to Gaussian
    model_output <- because_model(equations)

    expect_match(model_output$model, "dnorm")
    expect_false(grepl("dbern", model_output$model))
})

test_that("Poisson distribution with structure generates structural precision and sigma priors", {
    equations <- list(Abundance ~ X)
    distribution <- c(Abundance = "poisson")
    dummy_spatial <- diag(5)

    model_output <- because_model(
        equations,
        family = distribution,
        structures = list(spatial = dummy_spatial)
    )

    # Should define structural priors for Poisson response
    expect_match(model_output$model, "tau_u_spatial_Abundance")
    expect_match(model_output$model, "sigma_spatial_Abundance <- 1/sqrt\\(tau_u_spatial_Abundance\\)")
    # Should NOT declare residual variance tau_res for Poisson
    expect_false(grepl("tau_res_Abundance", model_output$model))
})
