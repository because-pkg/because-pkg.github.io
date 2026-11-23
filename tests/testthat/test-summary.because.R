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

test_that("summary.because prints cleanly without quotes when rope is specified", {
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
        dsep = FALSE
    )
    class(fit) <- "because"

    summ <- summary(fit, rope = c(-0.1, 0.1))
    output <- capture_output({
        print(summ)
    })

    expect_s3_class(summ, "summary.because")
    expect_true("P_ROPE" %in% colnames(summ$results))
    expect_false(grepl('"', output))
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
    expect_true(ncol(results$results) %in% c(8, 9))
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

test_that("summary.because filters redundant MAG indicator variance components", {
    # Mock samples with various variance parameters for MAG indicator variables
    param_names <- c(
        "beta_Brain_Mass",
        "rho_Brain_Mass",
        "sigma_res_Brain_Mass",
        "sigma_obs_Brain",
        "sigma_obs_Mass",
        "sigma_multiPhylo_Brain",
        "sigma_multiPhylo_Mass",
        "sigma_tot_Brain",
        "sigma_tot_Mass",
        "sigma_Brain_res",
        "sigma_Mass_res"
    )
    mat <- matrix(rnorm(100 * length(param_names), mean = 1, sd = 0.1), ncol = length(param_names), dimnames = list(NULL, param_names))
    samples <- coda::mcmc.list(coda::mcmc(mat))

    fit <- list(
        samples = samples,
        dsep = FALSE,
        induced_correlations = list(c("Brain", "Mass"))
    )
    class(fit) <- "because"

    summ <- summary(fit)
    rows <- rownames(summ$results)

    # Retained parameters
    expect_true("beta_Brain_Mass" %in% rows)
    expect_true("rho_Brain_Mass" %in% rows)
    expect_true("sigma_res_Brain_Mass" %in% rows)
    expect_true("sigma_obs_Brain" %in% rows)
    expect_true("sigma_obs_Mass" %in% rows)

    # Dropped redundant parameters
    expect_false("sigma_multiPhylo_Brain" %in% rows)
    expect_false("sigma_multiPhylo_Mass" %in% rows)
    expect_false("sigma_tot_Brain" %in% rows)
    expect_false("sigma_tot_Mass" %in% rows)
    expect_false("sigma_Brain_res" %in% rows)
    expect_false("sigma_Mass_res" %in% rows)
})

test_that("format_numpyro_samples retains >20 structural and grouped u_ vectors and notifies on omitted large vectors", {
    set.seed(42)
    py_res <- list(
        samples = list(
            alpha_Y = array(rnorm(2 * 10), dim = c(2, 10)),
            beta_Y_cat = array(rnorm(2 * 10 * 25), dim = c(2, 10, 25)),
            u_Y_Site = array(rnorm(2 * 10 * 35), dim = c(2, 10, 35)),
            u_std_Y_Site = array(rnorm(2 * 10 * 35), dim = c(2, 10, 35)),
            err_res_Y = array(rnorm(2 * 10 * 600), dim = c(2, 10, 600)),
            u_Y_phylo = array(rnorm(2 * 10 * 600), dim = c(2, 10, 600))
        )
    )

    # Under monitor = "interpretable", >20 structural and 35-group u_Y_Site are kept,
    # nuisance vectors are dropped, and omitting 600-species u_Y_phylo emits a clear message
    expect_message(
        mcmc_interp <- format_numpyro_samples(py_res, monitor = "interpretable", quiet = FALSE),
        "u_Y_phylo \\(length 600\\)"
    )
    cols_interp <- colnames(mcmc_interp[[1]])
    expect_true("alpha_Y" %in% cols_interp)
    expect_true("beta_Y_cat[25]" %in% cols_interp)
    expect_true("u_Y_Site[35]" %in% cols_interp)
    expect_false(any(grepl("^u_std_Y_Site\\[", cols_interp)))
    expect_false(any(grepl("^err_res_Y\\[", cols_interp)))
    expect_false(any(grepl("^u_Y_phylo\\[", cols_interp)))

    # Explicit monitor retains u_Y_phylo regardless of length without omission message
    expect_no_message(
        mcmc_custom <- format_numpyro_samples(py_res, monitor = c("interpretable", "u_Y_phylo"), quiet = FALSE)
    )
    expect_true("u_Y_phylo[600]" %in% colnames(mcmc_custom[[1]]))

    # monitor = "all" keeps everything including nuisance and high-dimensional vectors
    mcmc_all <- format_numpyro_samples(py_res, monitor = "all", quiet = FALSE)
    cols_all <- colnames(mcmc_all[[1]])
    expect_true("u_Y_phylo[600]" %in% cols_all)
    expect_true("u_std_Y_Site[35]" %in% cols_all)
    expect_true("err_res_Y[600]" %in% cols_all)
})

test_that("format_numpyro_samples handles 1-iteration and 1-element vector shapes safely", {
    py_edge <- list(
        samples = list(
            alpha_Y = array(1.5, dim = c(2, 1)),
            cutpoint_Y = array(0.25, dim = c(2, 1, 1)),
            u_Y_Site = array(c(0.1, -0.1, 0.2, -0.2), dim = c(2, 1, 2))
        )
    )
    mcmc_edge <- format_numpyro_samples(py_edge, monitor = "interpretable", quiet = TRUE)
    expect_s3_class(mcmc_edge, "mcmc.list")
    expect_equal(nrow(mcmc_edge[[1]]), 1)
    expect_true(all(c("alpha_Y", "cutpoint_Y[1]", "u_Y_Site[1]", "u_Y_Site[2]") %in% colnames(mcmc_edge[[1]])))
})

test_that("summary.because respects show_random = TRUE even with cached structural summary and warns when missing", {
    set.seed(123)
    mat <- cbind(
        alpha_Y = rnorm(40),
        beta_Y_X = rnorm(40, 0.5, 0.1),
        sigma_Y_Site = abs(rnorm(40, 0.4, 0.05)),
        `u_Y_Site[1]` = rnorm(40, -0.2, 0.1),
        `u_Y_Site[2]` = rnorm(40, 0.2, 0.1),
        `u_std_Y_Site[1]` = rnorm(40),
        `u_std_Y_Site[2]` = rnorm(40)
    )
    samples <- coda::mcmc.list(coda::mcmc(mat[1:20, ]), coda::mcmc(mat[21:40, ]))

    # Cached summary contains only structural parameters (simulating run_numpyro_pipeline)
    cached_sub <- samples[, c("alpha_Y", "beta_Y_X", "sigma_Y_Site")]
    cached_summ <- summary(cached_sub)
    cached_summ$statistics <- cbind(cached_summ$statistics, Rhat = c(1.0, 1.0, 1.0))

    fit <- list(
        samples = samples,
        summary = cached_summ,
        dsep = FALSE
    )
    class(fit) <- "because"

    # Default show_random = FALSE uses cached structural summary (no u_*)
    summ_def <- summary(fit)
    expect_false(any(grepl("^u_", rownames(summ_def$results))))

    # show_random = TRUE recomputes summary to expose u_Y_Site[*] while hiding auxiliary u_std_*
    summ_rand <- summary(fit, show_random = TRUE)
    expect_true("u_Y_Site[1]" %in% rownames(summ_rand$results))
    expect_true("u_Y_Site[2]" %in% rownames(summ_rand$results))
    expect_false("u_std_Y_Site[1]" %in% rownames(summ_rand$results))

    # Warning emitted when show_random = TRUE is requested on samples lacking u_*
    fit_no_u <- list(
        samples = cached_sub,
        summary = cached_summ,
        dsep = FALSE
    )
    class(fit_no_u) <- "because"
    expect_warning(
        summary(fit_no_u, show_random = TRUE),
        "No random effect estimates \\(u_\\*\\) were found"
    )
})
