test_that("Ordinal distribution runs successfully", {
    skip_on_cran()

    set.seed(456)
    N <- 100
    X <- rnorm(N)

    # Latent variable (Linear Predictor)
    # Note: eta = beta*X + err
    beta_val <- 1.5
    # Standard i.i.d. error instead of structured error
    err <- rnorm(N, 0, 0.5)
    eta <- beta_val * X + err

    # Cutpoints (K=3 categories)
    # P(Y <= 1) = logit^-1(c1 - eta)
    # P(Y <= 2) = logit^-1(c2 - eta)
    c1 <- -1.0
    c2 <- 1.0

    Q1 <- plogis(c1 - eta)
    Q2 <- plogis(c2 - eta)

    P1 <- Q1
    P2 <- Q2 - Q1
    P3 <- 1 - Q2

    Y <- numeric(N)
    for (i in 1:N) {
        Y[i] <- sample(1:3, 1, prob = c(P1[i], P2[i], P3[i]))
    }
    Y <- factor(Y, ordered = TRUE)

    data_list <- list(Y = Y, X = X)

    # Run model
    res <- because(
        equations = list(Y ~ X),
        data = data_list,
        
        family = c(Y = "ordinal"),
        n.iter = 1000,
        n.burnin = 500,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(res, "because")
    expect_true("beta_Y_X" %in% res$monitor)

    sum_stats <- res$summary$statistics

    # Check beta recovery
    # True value: 1.5
    expect_true("beta_Y_X" %in% rownames(sum_stats))
    beta_est <- sum_stats["beta_Y_X", "Mean"]
    expect_true(abs(beta_est - 1.5) < 1.0)

    # Check cutpoints
    # True values: -1, 1
    expect_true("cutpoint_Y[1]" %in% rownames(sum_stats))
    expect_true("cutpoint_Y[2]" %in% rownames(sum_stats))

    c1_est <- sum_stats["cutpoint_Y[1]", "Mean"]
    c2_est <- sum_stats["cutpoint_Y[2]", "Mean"]

    # Check relative order
    expect_true(c1_est < c2_est)
})
