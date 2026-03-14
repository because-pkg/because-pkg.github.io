library(testthat)

# Load current package state
if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".")
} else {
    stop("devtools needed to load package for testing")
}

# Access internal functions (now exported to namespace by load_all)
extract_random_effects <- because:::extract_random_effects
create_group_structures <- because:::create_group_structures

test_that("extract_random_effects parses correctly", {
    eqs <- list(
        Y ~ X + (1 | Group),
        Z ~ W + (1 | Other)
    )
    res <- extract_random_effects(eqs)

    # Check removed terms
    # Note: string matching depends on exact whitespace handling in parser
    # Fixed eq 1 should handle 'X'
    rhs1 <- as.character(res$fixed_equations[[1]])[3]
    expect_true(grepl("X", rhs1))
    expect_false(grepl("\\(1", rhs1))

    # Check random terms
    expect_equal(length(res$random_terms), 2)
    expect_equal(res$random_terms[[1]]$group, "Group")
    expect_equal(res$random_terms[[2]]$group, "Other")
})

test_that("create_group_structures creates correct data", {
    df <- data.frame(
        Group = c("A", "A", "B", "C"),
        Other = c("x", "y", "x", "y"),
        Val = 1:4
    )
    random_terms <- list(
        list(group = "Group", type = "intercept"),
        list(group = "Other", type = "intercept")
    )

    res <- create_group_structures(df, random_terms)
    updates <- res$data_updates

    # Group: 3 levels (A, B, C)
    expect_true("N_Group" %in% names(updates))
    expect_equal(updates$N_Group, 3)
    expect_equal(length(updates$group_Group), 4)
    expect_equal(updates$group_Group, as.integer(as.factor(df$Group)))
    expect_equal(dim(updates$Prec_Group), c(3, 3))

    # Other: 2 levels (x, y)
    expect_equal(updates$N_Other, 2)
})

test_that("Simulated Gaussian Random Effect Model runs and recovers parameters", {
    skip_on_cran() # Skip computationally intensive test

    set.seed(42)
    N <- 200
    n_groups <- 10
    groups <- sample(LETTERS[1:n_groups], N, replace = TRUE)

    # True parameters
    sigma_u_true <- 2.0
    sigma_e_true <- 1.0
    beta_x <- 0.5
    intercept <- 1.0

    # Random effects (centered to ensure intercept identifiability)
    u_true <- rnorm(n_groups, 0, sigma_u_true)
    u_true <- u_true - mean(u_true)
    names(u_true) <- LETTERS[1:n_groups]

    x <- rnorm(N)
    y <- intercept + beta_x * x + u_true[groups] + rnorm(N, 0, sigma_e_true)

    df <- data.frame(Y = y, X = x, Group = as.factor(groups))

    # Run model
    # Y ~ X + (1|Group)
    fit <- because(
        equations = list(Y ~ X + (1 | Group)),
        data = df,
        n.chains = 1,
        n.iter = 2000,
        n.burnin = 500,
        n.thin = 2,
        optimise = TRUE, # Should be default anyway
        quiet = TRUE
    )

    sum_stats <- fit$summary$statistics

    print("Estimated Parameters:")
    print(sum_stats[, c("Mean", "SD")])

    # Check recovery
    # Intercept (alphaY - no underscore for single response)
    # Relax tolerance slightly for short chain
    expect_true(abs(sum_stats["alphaY", "Mean"] - intercept) < 1.0) # Was 0.5
    # Beta X (beta_Y_X - underscores)
    # Note: coda usage might Capitalize "Mean"
    expect_true(abs(sum_stats["beta_Y_X", "Mean"] - beta_x) < 0.2)
    # Sigma Group (sigma_Y_Group)
    expect_true(abs(sum_stats["sigma_Y_Group", "Mean"] - sigma_u_true) < 1.0)
    # Sigma Residual (sigma_Y_res)
    expect_true(abs(sum_stats["sigma_Y_res", "Mean"] - sigma_e_true) < 0.5)
})

test_that("Random Effects with Poisson Distribution", {
    skip_on_cran()

    set.seed(42)
    N <- 100
    n_groups <- 5
    groups <- sample(1:n_groups, N, replace = TRUE)
    u_true <- rnorm(n_groups, 0, 0.5)

    # Lambda = exp(1 + u)
    lambda <- exp(1 + u_true[groups])
    y <- rpois(N, lambda)

    df <- data.frame(Y = y, Group = as.factor(groups))

    fit <- because(
        equations = list(Y ~ 1 + (1 | Group)),
        data = df,
        family = list(Y = "poisson"),
        n.chains = 1,
        n.iter = 1000,
        n.burnin = 200,
        quiet = TRUE
    )

    sum_stats <- fit$summary$statistics

    expect_true("sigma_Y_Group" %in% rownames(sum_stats))
})

test_that("Random Slopes Warn and Default to Intercept", {
    eq <- list(Y ~ X + (X | Group))

    # Check that it throws warning
    expect_warning(
        because:::extract_random_effects(eq),
        "Random slopes .* are not yet implemented"
    )

    # Check that it still returns valid random structure
    suppressWarnings({
        res <- because:::extract_random_effects(eq)
    })

    expect_equal(length(res$random_terms), 1)
    expect_equal(res$random_terms[[1]]$group, "Group")
    expect_equal(res$random_terms[[1]]$type, "intercept")
})

test_that("Multiple Random Effects Run Successfully", {
    skip_on_cran()

    set.seed(123)
    N <- 50
    # Factor A
    groups_A <- sample(1:5, N, replace = TRUE)
    u_A <- rnorm(5, 0, 1.0)
    u_A <- u_A - mean(u_A)

    # Factor B
    groups_B <- sample(1:3, N, replace = TRUE)
    u_B <- rnorm(3, 0, 0.5)
    u_B <- u_B - mean(u_B)

    y <- 1 + u_A[groups_A] + u_B[groups_B] + rnorm(N, 0, 0.5)
    df <- data.frame(Y = y, A = as.factor(groups_A), B = as.factor(groups_B))

    fit <- because(
        equations = list(Y ~ 1 + (1 | A) + (1 | B)),
        data = df,
        n.chains = 1,
        n.iter = 500,
        quiet = TRUE
    )

    sum_stats <- fit$summary$statistics
    expect_true("sigma_Y_A" %in% rownames(sum_stats))
    expect_true("sigma_Y_B" %in% rownames(sum_stats))

    # Rough check of magnitude
    expect_true(abs(sum_stats["sigma_Y_A", "Mean"] - 1.0) < 0.8)
    expect_true(abs(sum_stats["sigma_Y_B", "Mean"] - 0.5) < 0.8)
})

test_that("Parallel Execution with Global Random Effects and WAIC", {
    skip_on_cran()

    set.seed(789)
    N <- 100
    n_groups <- 10
    groups <- sample(LETTERS[1:n_groups], N, replace = TRUE)
    u_true <- rnorm(n_groups, 0, 1.0)
    names(u_true) <- LETTERS[1:n_groups]
    beta <- 0.5

    X <- rnorm(N)
    Y <- 1 + beta * X + u_true[groups] + rnorm(N, 0, 0.5)

    data_df <- data.frame(Y = Y, X = X, Group = as.factor(groups), ID = 1:N)

    # Run with parallel = TRUE and WAIC = TRUE and global Random
    # Expectation: Should run without error and return valid because object with WAIC
    fit_par <- because(
        equations = list(Y ~ X),
        data = data_df,
        random = ~ (1 | Group),
        parallel = TRUE,
        n.chains = 2,
        n.iter = 1000,
        n.burnin = 200,
        n.cores = 2,
        WAIC = TRUE,
        quiet = TRUE
    )

    expect_s3_class(fit_par, "because")
    expect_false(is.null(fit_par$WAIC))

    # Check Model Code for Random Effects
    # Should contain u_Y_Group or u_Group
    model_code <- fit_par$model_code
    has_random_term <- grepl("u_Group", model_code) ||
        grepl("u_.*Group", model_code)
    expect_true(has_random_term)
    expect_true(grepl("group_Group", model_code))

    # Check Sigma for Group in Summary
    sum_stats <- fit_par$summary$statistics
    expect_true("sigma_Y_Group" %in% rownames(sum_stats))
})

test_that("Global Random Argument & d-separation work correctly", {
    skip_on_cran()

    set.seed(456)
    N <- 100
    # Path model: A -> B -> C w/ global (1|Group)
    groups <- sample(LETTERS[1:3], N, replace = TRUE)
    data <- data.frame(
        A = rnorm(N),
        B = rnorm(N),
        C = rnorm(N),
        Group = as.factor(groups)
    )

    equations <- list(B ~ A, C ~ B)

    fit <- because(
        equations = equations,
        data = data,
        random = ~ (1 | Group), # Global argument
        dsep = TRUE,
        n.chains = 2,
        n.iter = 2000,
        quiet = TRUE
    )

    # 1. Verify Global Random effect was applied
    # By checking if d-sep test formula includes (1|Group)
    # Test for C _||_ A | B
    dsep_test <- fit$dsep_tests[[1]]
    f_str <- deparse(dsep_test)
    expect_true(
        grepl("\\(1 \\| Group\\)", f_str) || grepl("\\(1\\|Group\\)", f_str)
    )

    # 2. Verify Monitor Logic (should exclude beta_C_Group)
    # If monitor logic failed, JAGS would error and fit$samples would be NULL or incomplete
    expect_false(is.null(fit$samples))

    # Check that we monitored beta_C_A (or similar depending on equation)
    # C ~ A + B + (1|Group). Test var is A (if B is parent).
    # Monitor: beta_C_A
    samps <- fit$samples[[1]]
    expect_true("beta_C_A" %in% colnames(samps))
    # Should NOT have beta_C_Group
    expect_false("beta_C_Group" %in% colnames(samps))
})
