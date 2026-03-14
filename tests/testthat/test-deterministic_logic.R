context("Deterministic Logic Operators")

test_that("term_to_jags_expression converts R logic to JAGS logic", {
    # We access the internal function
    # Note: running this requires 'because' to be loaded.
    # If running via devtools::test(), it is loaded.

    f <- because:::term_to_jags_expression

    # Basic variable expansion + logic replacement
    expect_equal(f("A & B"), "A[i] && B[i]")
    expect_equal(f("A | B"), "A[i] || B[i]")

    # Preservation of other operators
    expect_equal(f("A >= 1 & A < 2"), "A[i] >= 1 && A[i] < 2")
    expect_equal(f("(A > 0) & (B < 0)"), "(A[i] > 0) && (B[i] < 0)")

    # Multiple operators
    expect_equal(f("A & B & C"), "A[i] && B[i] && C[i]")

    # Interaction with == replacement
    # A == 1 & B == 2 -> equals(A[i], 1) && equals(B[i], 2)
    expect_equal(f("A == 1 & B == 2"), "equals(A[i], 1) && equals(B[i], 2)")
})

test_that("Model fits with logical AND in formula", {
    skip_on_cran()

    set.seed(123)
    N <- 50
    A <- runif(N, 0, 10)
    # Define Y based on range of A
    # logic: if A in [2, 5], Y increases by 2
    in_range <- (A >= 2 & A < 5)
    Y <- A + 2 * in_range + rnorm(N, sd = 0.1)

    data <- list(Y = Y, A = A)

    # Formula using &
    eqs <- list(
        Y ~ A + I(A >= 2 & A < 5)
    )

    # We expect this to run without syntax error in JAGS
    # Suppress output for clean logs
    fit <- because(
        equations = eqs,
        data = data,
        n.iter = 50,
        quiet = TRUE,
        n.chains = 1
    )

    # Check if the interaction term exists in results
    stats <- fit$summary$statistics

    # Identify coefficients
    coefs <- rownames(stats)
    # Expect beta_Y_A_gt_2_and_A_lt_5
    expect_true(any(grepl("beta_Y_A_gt_2", coefs)))

    # Check estimate (should be approx 2)
    # finding the specific one
    det_coef <- coefs[grepl("beta_Y_A_gt_2", coefs)]
    val <- stats[det_coef, "Mean"]
    expect_equal(val, 2, tolerance = 0.5)
})
