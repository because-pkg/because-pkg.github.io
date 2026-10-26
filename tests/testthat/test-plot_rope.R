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
