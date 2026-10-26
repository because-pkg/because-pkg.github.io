test_that("Multinomial distribution runs successfully", {
    skip_on_cran()

    set.seed(123)
    N <- 100
    K <- 3
    X <- rnorm(N)

    # Parameters (Cat 1 reference)
    # Latent 1 (Cat 2 vs 1)
    beta1_0 <- 0.5
    beta1_1 <- 1.2
    phy1 <- rnorm(N, sd = 0.5)
    L1 <- beta1_0 + beta1_1 * X + phy1

    # Latent 2 (Cat 3 vs 1)
    beta2_0 <- -0.5
    beta2_1 <- -0.8
    phy2 <- rnorm(N, sd = 0.5)
    L2 <- beta2_0 + beta2_1 * X + phy2

    # Softmax
    P <- matrix(0, N, K)
    denom <- 1 + exp(L1) + exp(L2)
    P[, 1] <- 1 / denom
    P[, 2] <- exp(L1) / denom
    P[, 3] <- exp(L2) / denom

    Y <- numeric(N)
    for (i in 1:N) {
        Y[i] <- sample(1:K, 1, prob = P[i, ])
    }
    Y <- factor(Y) # Ensure it's a factor

    data_list <- list(Y = Y, X = X)

    # Run model
    res <- because(
        equations = list(Y ~ X),
        data = data_list,
        
        family = c(Y = "multinomial"),
        n.iter = 1000,
        n.burnin = 500,
        n.chains = 2,
        quiet = TRUE
    )

    expect_s3_class(res, "because")
    expect_true("beta_Y_X" %in% res$monitor)

    # Check summary
    sum_stats <- res$summary$statistics

    # We expect beta_Y_X[2] and beta_Y_X[3]
    # Note: JAGS indexing might be beta_Y_X[2] and beta_Y_X[3] because loop is 2:K
    # Let's check the row names
    expect_true(any(grepl("beta_Y_X\\[2\\]", rownames(sum_stats))))
    expect_true(any(grepl("beta_Y_X\\[3\\]", rownames(sum_stats))))

    # Check recovery (loose bounds due to small N)
    b1_est <- sum_stats["beta_Y_X[2]", "Mean"]
    b2_est <- sum_stats["beta_Y_X[3]", "Mean"]

    # True values: 1.2 and -0.8
    expect_true(abs(b1_est - 1.2) < 1.5) # Loose check
    expect_true(abs(b2_est - (-0.8)) < 1.5)
})
