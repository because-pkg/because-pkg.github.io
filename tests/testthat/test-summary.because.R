test_that("summary.because handles standard output", {
    # Mock a because object
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(100),
            ncol = 2,
            dimnames = list(NULL, c("alpha", "beta"))
        )),
        coda::mcmc(matrix(
            rnorm(100),
            ncol = 2,
            dimnames = list(NULL, c("alpha", "beta"))
        ))
    )

    fit <- list(
        samples = samples,
        dsep = FALSE,
        DIC = c(deviance = 100, pD = 2, DIC = 102),
        WAIC = c(waic = 100, p_waic = 2)
    )
    class(fit) <- "because"

    # Capture output to verify printing (explicit print required)
    summ <- summary(fit)
    output <- capture_output({
        print(summ)
    })

    # summary.because returns a summary.because object
    expect_s3_class(summ, "summary.because")
    expect_true(is.matrix(summ$results) || is.data.frame(summ$results))
    expect_match(output, "DIC")
    expect_match(output, "WAIC")
})

test_that("summary.because handles d-sep output", {
    # Mock a d-sep because object
    # Use 2 parameters to ensure summary returns a matrix
    samples <- coda::mcmc.list(
        coda::mcmc(matrix(
            rnorm(400, mean = 0),
            ncol = 2,
            dimnames = list(NULL, c("beta_Y_X", "beta_Z_Y"))
        ))
    )

    # Setup d-sep specific fields
    dsep_tests <- list(Y ~ X)
    attr(dsep_tests[[1]], "test_var") <- "X"

    parameter_map <- data.frame(
        response = "Y",
        predictor = "X",
        parameter = "beta_Y_X",
        equation_index = 1,
        stringsAsFactors = FALSE
    )

    fit <- list(
        samples = samples,
        dsep = TRUE,
        dsep_tests = dsep_tests,
        parameter_map = parameter_map,
        dsep_results = list(
            list(
                samples = samples,
                param_map = parameter_map
            )
        )
    )
    class(fit) <- "because"

    results <- summary(fit)
    output <- capture_output({
        print(results)
    })

    # Output matches "d-separation Tests" (updated string)
    expect_match(output, "d-separation Tests")
    expect_match(output, "Y _\\|\\|_ X")

    # Check structure of returned object
    expect_s3_class(results, "summary.because")
    expect_true(is.data.frame(results$results))
    expect_true(nrow(results$results) == 1)
    expect_equal(ncol(results$results), 8)
    expect_true("LowerCI" %in% colnames(results$results))
})

test_that("format_numpyro_samples extracts rho from 4D L_corr and summary retains rho", {
    # Mock NumPyro output with 2 chains, 50 iterations
    # 2D scalar: beta_y_x
    beta_data <- array(rnorm(100), dim = c(2, 50))
    # 4D Cholesky factor: L_corr_Brain_Mass (2 chains, 50 iters, 2x2)
    L_corr_data <- array(0, dim = c(2, 50, 2, 2))
    # Diagonal element (1,1) is 1.0, (2,2) is sqrt(1 - rho^2)
    # Off-diagonal (2,1) is rho
    L_corr_data[, , 1, 1] <- 1.0
    L_corr_data[, , 2, 1] <- 0.75
    L_corr_data[, , 2, 2] <- sqrt(1 - 0.75^2)

    py_res <- list(
        samples = list(
            beta_y_x = beta_data,
            L_corr_Brain_Mass = L_corr_data
        )
    )

    mcmc_list <- format_numpyro_samples(py_res)
    expect_s3_class(mcmc_list, "mcmc.list")
    expect_equal(length(mcmc_list), 2)
    expect_true("rho_Brain_Mass" %in% colnames(mcmc_list[[1]]))
    expect_equal(unname(mcmc_list[[1]][1, "rho_Brain_Mass"]), 0.75)

    fit <- list(
        samples = mcmc_list,
        dsep = FALSE
    )
    class(fit) <- "because"

    summ <- summary(fit)
    expect_true("rho_Brain_Mass" %in% rownames(summ$results))
    expect_equal(round(summ$results["rho_Brain_Mass", "Mean"], 2), 0.75)
})
