library(because)

test_that("standardize works on fitted because models with Gaussian and Binomial traits", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    y_raw <- rnorm(n, mean = 2 * x, sd = 2) # SD(Y) ~ 2.8, SD(X) ~ 1
    migrant <- rbinom(n, 1, 0.4) # binary 0/1 with SD ~ sqrt(0.24) ~ 0.49

    data_df <- data.frame(
        X = x,
        Y = y_raw,
        Migrant = migrant,
        Species = paste0("sp_", 1:n)
    )

    # Mock samples for 2 chains, 200 iterations
    # beta_Y_X ~ 2.0 (unstandardized) -> standardized beta should be 2.0 * (1 / sd(Y)) ~ 0.7
    # beta_Migrant_X ~ 1.5 (logit scale)
    # sigma_Migrant_phylo ~ 1.2
    # beta_Y_Migrant ~ 1.0 (difference between migrant 1 and 0) -> standardized should be 1.0 * (sd(Migrant) / sd(Y))
    n_iter <- 200
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            c(
                rnorm(n_iter, 2.0, 0.05), # beta_Y_X
                rnorm(n_iter, 1.0, 0.05), # beta_Y_Migrant
                rnorm(n_iter, 1.5, 0.05), # beta_Migrant_X
                abs(rnorm(n_iter, 1.2, 0.05)) # sigma_Migrant_phylo
            ),
            ncol = 4,
            dimnames = list(NULL, c("beta_Y_X", "beta_Y_Migrant", "beta_Migrant_X", "sigma_Migrant_phylo"))
        )),
        coda::mcmc(matrix(
            c(
                rnorm(n_iter, 2.0, 0.05),
                rnorm(n_iter, 1.0, 0.05),
                rnorm(n_iter, 1.5, 0.05),
                abs(rnorm(n_iter, 1.2, 0.05))
            ),
            ncol = 4,
            dimnames = list(NULL, c("beta_Y_X", "beta_Y_Migrant", "beta_Migrant_X", "sigma_Migrant_phylo"))
        ))
    )

    param_map <- data.frame(
        response = c("Y", "Y", "Migrant"),
        predictor = c("X", "Migrant", "X"),
        parameter = c("beta_Y_X", "beta_Y_Migrant", "beta_Migrant_X"),
        equation_index = c(1, 1, 2),
        stringsAsFactors = FALSE
    )

    fit <- list(
        samples = samples,
        parameter_map = param_map,
        data = data_df,
        family = c(Y = "gaussian", Migrant = "binomial"),
        dsep = FALSE,
        equations = list(Y ~ X + Migrant, Migrant ~ X)
    )
    class(fit) <- "because"

    # Standardize
    fit_std <- standardize(fit)

    expect_s3_class(fit_std, "because_standardized")
    expect_s3_class(fit_std, "because")
    expect_true(isTRUE(attr(fit_std, "standardized")))

    # Verify Gaussian rescaling: beta_Y_X was multiplied by sd(X)/sd(Y)
    sd_x <- sd(data_df$X)
    sd_y <- sd(data_df$Y)
    expected_beta_Y_X <- mean(samples[[1]][, "beta_Y_X"]) * (sd_x / sd_y)
    expect_equal(mean(fit_std$samples[[1]][, "beta_Y_X"]), expected_beta_Y_X, tolerance = 1e-4)

    # Verify binary predictor into Gaussian: beta_Y_Migrant was multiplied by sd(Migrant)/sd(Y)
    sd_mig <- sd(data_df$Migrant)
    expected_beta_Y_Mig <- mean(samples[[1]][, "beta_Y_Migrant"]) * (sd_mig / sd_y)
    expect_equal(mean(fit_std$samples[[1]][, "beta_Y_Migrant"]), expected_beta_Y_Mig, tolerance = 1e-4)

    # Verify Binomial latent liability rescaling: beta_Migrant_X was divided by sqrt(sigma_phylo^2 + pi^2/3)
    sig_phylo <- samples[[1]][, "sigma_Migrant_phylo"]
    sd_liability <- sqrt(sig_phylo^2 + (pi^2) / 3)
    expected_beta_Mig_X <- mean(samples[[1]][, "beta_Migrant_X"] * (sd_x / sd_liability))
    expect_equal(mean(fit_std$samples[[1]][, "beta_Migrant_X"]), expected_beta_Mig_X, tolerance = 1e-4)

    # Idempotence: calling standardize again does not double-rescale
    fit_std2 <- standardize(fit_std)
    expect_equal(mean(fit_std2$samples[[1]][, "beta_Y_X"]), mean(fit_std$samples[[1]][, "beta_Y_X"]))

    # Summary works and displays note
    s <- summary(fit_std)
    expect_true(isTRUE(s$standardized))
    out <- capture.output(print(s))
    expect_true(any(grepl("Model is standardized", out)))
})

test_that("standardize works on d-separation objects and harmonizes ROPE", {
    set.seed(42)
    n <- 100
    x <- rnorm(n)
    migrant <- rbinom(n, 1, 0.4)
    y <- rnorm(n)

    data_df <- data.frame(
        X = x,
        Migrant = migrant,
        Y = y
    )

    n_iter <- 150
    # Test 1: Y _||_ Migrant | {}  (Y Gaussian, Migrant binary)
    # beta_Y_Migrant ~ 0.5
    samples1 <- coda::mcmc.list(
        coda::mcmc(matrix(rnorm(n_iter, 0.5, 0.05), ncol = 1, dimnames = list(NULL, "beta_Y_Migrant"))),
        coda::mcmc(matrix(rnorm(n_iter, 0.5, 0.05), ncol = 1, dimnames = list(NULL, "beta_Y_Migrant")))
    )
    test1_formula <- Y ~ Migrant
    attr(test1_formula, "test_var") <- "Migrant"

    # Test 2: Migrant _||_ X | {} (Migrant Binomial, X continuous)
    # beta_Migrant_X ~ 0.8
    samples2 <- coda::mcmc.list(
        coda::mcmc(matrix(rnorm(n_iter, 0.8, 0.05), ncol = 1, dimnames = list(NULL, "beta_Migrant_X"))),
        coda::mcmc(matrix(rnorm(n_iter, 0.8, 0.05), ncol = 1, dimnames = list(NULL, "beta_Migrant_X")))
    )
    test2_formula <- Migrant ~ X
    attr(test2_formula, "test_var") <- "X"

    combined_samples <- coda::mcmc.list(
        coda::mcmc(cbind(as.matrix(samples1[[1]]), as.matrix(samples2[[1]]))),
        coda::mcmc(cbind(as.matrix(samples1[[2]]), as.matrix(samples2[[2]])))
    )

    dsep_fit <- list(
        samples = combined_samples,
        dsep = TRUE,
        dsep_tests = list(test1_formula, test2_formula),
        dsep_results = list(
            list(samples = samples1, param_map = data.frame(response = "Y", predictor = "Migrant", parameter = "beta_Y_Migrant", equation_index = 1)),
            list(samples = samples2, param_map = data.frame(response = "Migrant", predictor = "X", parameter = "beta_Migrant_X", equation_index = 1))
        ),
        data = data_df,
        family = c(Y = "gaussian", Migrant = "binomial")
    )
    class(dsep_fit) <- "because"

    # Standardize
    dsep_std <- standardize(dsep_fit)

    expect_s3_class(dsep_std, "because_standardized")
    expect_true(isTRUE(attr(dsep_std, "standardized")))

    # In test 1: beta_Y_Migrant is scaled by sd(Migrant)/sd(Y)
    sd_m <- sd(data_df$Migrant)
    sd_y <- sd(data_df$Y)
    expected_b1 <- mean(samples1[[1]][, "beta_Y_Migrant"]) * (sd_m / sd_y)
    expect_equal(mean(dsep_std$dsep_results[[1]]$samples[[1]][, "beta_Y_Migrant"]), expected_b1, tolerance = 1e-4)

    # In test 2: beta_Migrant_X is scaled by sd(X) / (pi / sqrt(3))
    sd_x <- sd(data_df$X)
    expected_b2 <- mean(samples2[[1]][, "beta_Migrant_X"]) * (sd_x / (pi / sqrt(3)))
    expect_equal(mean(dsep_std$dsep_results[[2]]$samples[[1]][, "beta_Migrant_X"]), expected_b2, tolerance = 1e-4)

    # Global samples are updated as well
    expect_equal(mean(dsep_std$samples[[1]][, "beta_Y_Migrant"]), mean(dsep_std$dsep_results[[1]]$samples[[1]][, "beta_Y_Migrant"]))
    expect_equal(mean(dsep_std$samples[[1]][, "beta_Migrant_X"]), mean(dsep_std$dsep_results[[2]]$samples[[1]][, "beta_Migrant_X"]))

    # Summary with rope = TRUE uses c(-0.10, 0.10) uniformly across both tests
    s <- summary(dsep_std, rope = TRUE)
    expect_true(isTRUE(s$standardized))
    out <- capture.output(print(s))
    expect_true(any(grepl("Note: d-separation tests are evaluated on the standardized scale", out)))

    # plot_rope works on standardized dsep object
    p_rope <- plot_rope(dsep_std)
    expect_s3_class(p_rope, "ggplot")

    # Both tests in plot_rope have r_min = -0.1 and r_max = 0.1
    expect_equal(unique(p_rope$data$r_min), -0.1)
    expect_equal(unique(p_rope$data$r_max), 0.1)
})

test_that("plot_coef and plot_dag work with standardized objects", {
    set.seed(123)
    n <- 50
    x <- rnorm(n)
    y <- rnorm(n)
    data_df <- data.frame(X = x, Y = y, Species = paste0("sp", 1:n))
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(rnorm(100, 0.4, 0.05), ncol = 1, dimnames = list(NULL, "beta_Y_X"))),
        coda::mcmc(matrix(rnorm(100, 0.4, 0.05), ncol = 1, dimnames = list(NULL, "beta_Y_X")))
    )
    param_map <- data.frame(
        response = "Y",
        predictor = "X",
        parameter = "beta_Y_X",
        equation_index = 1,
        stringsAsFactors = FALSE
    )
    fit <- list(
        samples = samples,
        parameter_map = param_map,
        data = data_df,
        family = c(Y = "gaussian"),
        equations = list(Y ~ X)
    )
    class(fit) <- "because"

    fit_std <- standardize(fit)

    # plot_coef works
    p_coef <- plot_coef(fit_std)
    expect_s3_class(p_coef, "ggplot")

    # plot_dag works
    p_dag <- plot_dag(fit_std)
    expect_s3_class(p_dag, "ggplot")
})

test_that("standardize correctly detects binomial from 0/1 data when family attribute is omitted", {
    set.seed(123)
    n <- 50
    x <- rnorm(n)
    migrant <- rbinom(n, 1, 0.5) # binary 0/1

    data_df <- data.frame(
        Brain = x,
        Migrant = migrant,
        Species = paste0("sp_", 1:n)
    )

    n_iter <- 100
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            c(
                rnorm(n_iter, -1.4, 0.05), # beta_Migrant_Brain
                rep(3.4, n_iter)           # sigma_Migrant_multiPhylo
            ),
            ncol = 2,
            dimnames = list(NULL, c("beta_Migrant_Brain", "sigma_Migrant_multiPhylo"))
        ))
    )

    param_map <- data.frame(
        response = c("Migrant"),
        predictor = c("Brain"),
        parameter = c("beta_Migrant_Brain"),
        equation_index = c(1),
        stringsAsFactors = FALSE
    )

    # Note: family is intentionally NULL to test fallback detection
    fit <- list(
        samples = samples,
        parameter_map = param_map,
        data = data_df,
        family = NULL,
        equations = list(Migrant ~ Brain)
    )
    class(fit) <- "because"

    fit_std <- standardize(fit)

    # Expected liability SD: sqrt(3.4^2 + pi^2/3) ~ sqrt(11.56 + 3.289868) ~ 3.85355
    expected_sd_liab <- sqrt(3.4^2 + (pi^2) / 3)
    sd_brain <- sd(data_df$Brain)
    expected_beta <- -1.4 * (sd_brain / expected_sd_liab)

    # Must be around -0.36, definitely NOT doubled (-2.8)
    actual_beta <- mean(fit_std$samples[[1]][, "beta_Migrant_Brain"])
    expect_equal(actual_beta, expected_beta, tolerance = 0.05)
    expect_true(actual_beta > -0.5 && actual_beta < -0.2)
})

