test_that("plot_rope filters to betas and rhos by default on fitted models", {
    # Mock samples containing betas, alphas, and sigmas
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(600, mean = 0.5),
            ncol = 6,
            dimnames = list(NULL, c(
                "beta_prey_habitat", "beta_abundance_prey",
                "alpha_prey", "alpha_abundance",
                "sigma_prey_res", "sigma_abundance_res"
            ))
        ))
    )

    parameter_map <- data.frame(
        response = c("prey", "abundance"),
        predictor = c("habitat", "prey"),
        parameter = c("beta_prey_habitat", "beta_abundance_prey"),
        equation_index = c(1, 2),
        stringsAsFactors = FALSE
    )

    summ <- summary(samples)

    fit <- list(
        samples = samples,
        dsep = FALSE,
        parameter_map = parameter_map,
        summary = summ
    )
    class(fit) <- "because"

    # Only run plot test if ggridges is available
    skip_if_not_installed("ggridges")

    p <- plot_rope(fit)

    expect_s3_class(p, "ggplot")
    # By default, only the two betas should be plotted with friendly labels
    expect_equal(levels(p$data$Label), c("abundance ~ prey", "prey ~ habitat"))

    # Explicit override: parameters = "all" includes alpha and sigma
    p_all <- plot_rope(fit, parameters = "all")
    expect_equal(length(unique(p_all$data$Label)), 6)
})

test_that("plot_rope disambiguates categorical dummy variables for d-sep tests", {
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(600, mean = 0),
            ncol = 2,
            dimnames = list(NULL, c("beta_Brain_Sample_size_medium", "beta_Brain_Sample_size_small"))
        ))
    )

    dsep_tests <- list(Brain ~ Sample.size)
    attr(dsep_tests[[1]], "test_var") <- "Sample.size"
    attr(dsep_tests[[1]], "formula_string") <- "Brain _||_ Sample.size"

    param_map <- data.frame(
        response = c("Brain", "Brain"),
        predictor = c("Sample_size_medium", "Sample_size_small"),
        parameter = c("beta_Brain_Sample_size_medium", "beta_Brain_Sample_size_small"),
        equation_index = c(1, 1),
        stringsAsFactors = FALSE
    )

    fit <- list(
        samples = samples,
        dsep = TRUE,
        dsep_tests = dsep_tests,
        parameter_map = param_map,
        dsep_results = list(
            list(
                samples = samples,
                param_map = param_map
            )
        )
    )
    class(fit) <- "because"

    skip_if_not_installed("ggridges")
    p <- plot_rope(fit)

    expect_s3_class(p, "ggplot")
    labels <- levels(p$data$Label)
    expect_true("Brain _||_ Sample.size (medium)" %in% labels)
    expect_true("Brain _||_ Sample.size (small)" %in% labels)
})
