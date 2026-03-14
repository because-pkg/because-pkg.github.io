context("Deterministic Nodes")

test_that("Deterministic Nodes handle Interactions, Logic, and Math", {
    skip_on_cran()

    set.seed(123)
    N <- 50
    tree <- ape::rtree(N)

    # Predictors
    A <- runif(N, min = 1, max = 5) # Ensure strictly positive for log()
    B <- rnorm(N)

    # 1. Interaction: Y_int ~ A * B
    Y_int <- 1 + 2 * (A * B) + rnorm(N, sd = 0.1)

    # 2. Logic: Y_logic ~ (A > 2)
    # JAGS (A > 2) returns 1 if true, 0 if false
    Y_logic <- 1 + 3 * as.numeric(A > 2) + rnorm(N, sd = 0.1)

    # 3. Math: Y_math ~ log(A)
    Y_math <- 1 + 4 * log(A) + rnorm(N, sd = 0.1)

    # 4. Arithmetic: Y_sum ~ (A + B)
    Y_sum <- 1 + 5 * (A + B) + rnorm(N, sd = 0.1)

    data <- list(
        Y_int = Y_int,
        Y_logic = Y_logic,
        Y_math = Y_math,
        Y_sum = Y_sum,
        A = A,
        B = B
    )

    # Fit Model
    # We test all distinct types in one model to save time
    # Note: logic requires I() wrapper in formula
    equations <- list(
        Y_int ~ A * B,
        Y_logic ~ I(A > 2),
        Y_math ~ I(log(A)),
        Y_sum ~ I(A + B)
    )

    fit <- because(
        equations = equations,
        data = data,

        n.iter = 100,
        n.chains = 1,
        quiet = FALSE
    )

    print("DEBUG: Checking fit object structure")
    print(str(fit))
    stats <- fit$summary$statistics

    # Verification Logic

    # 1. Interaction
    # Internal name should be A_x_B
    # Coefficient: beta_Y_int_A_x_B
    beta_int <- stats[rownames(stats) == "beta_Y_int_A_x_B", "Mean"]
    expect_true(!is.na(beta_int))
    expect_equal(beta_int, 2, tolerance = 0.5)

    # 2. Logic
    # Internal name: A_gt_2 (sanitized A > 2)
    # Coefficient: beta_Y_logic_A_gt_2
    param_name <- grep("beta_Y_logic_A_gt_2", rownames(stats), value = TRUE)
    beta_logic <- stats[param_name, "Mean"]
    expect_true(!is.na(beta_logic))
    expect_equal(beta_logic, 3, tolerance = 0.5)

    # 3. Math
    # Internal name: log_A_ (sanitized log(A)) - parentheses removed/replaced
    # Let's check sanitization: log(A) -> log_div_A_div_ -> log_A
    # wait, sanitize: () not handled explicitly in my function except generally?
    # sanitize: gsub("[^a-zA-Z0-9_]", "_", out)
    # log(A) -> log_A
    beta_math <- stats[rownames(stats) == "beta_Y_math_log_A", "Mean"]
    expect_true(!is.na(beta_math))
    expect_equal(beta_math, 4, tolerance = 0.5)

    # I(A+B) -> A_plus_B
    param_name_sum <- grep("beta_Y_sum_A_plus_B", rownames(stats), value = TRUE)
    beta_sum <- stats[param_name_sum, "Mean"]
    expect_true(!is.na(beta_sum))
    expect_equal(beta_sum, 5, tolerance = 0.5)
})
