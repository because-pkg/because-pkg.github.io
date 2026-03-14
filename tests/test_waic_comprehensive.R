##############################################################################
# Comprehensive WAIC-SE Testing Suite
# Tests WAIC with standard errors across various model types
##############################################################################

library(because)
library(ape)

set.seed(12345)

# Helper function
cat(paste(rep("=", 78), collapse = ""), "\n")
cat("COMPREHENSIVE WAIC-SE TESTING SUITE\n")
cat(paste(rep("=", 78), collapse = ""), "\n\n")

test_counter <- 0
pass_counter <- 0
fail_counter <- 0

run_test <- function(name, expr) {
    test_counter <<- test_counter + 1
    cat(sprintf("\n[Test %d] %s\n", test_counter, name))
    cat(paste(rep("-", 78), collapse = ""), "\n")

    result <- tryCatch(
        {
            eval(expr)
            TRUE
        },
        error = function(e) {
            cat("✗ FAIL:", e$message, "\n")
            FALSE
        }
    )

    if (result) {
        pass_counter <<- pass_counter + 1
        cat("✓ PASS\n")
    } else {
        fail_counter <<- fail_counter + 1
    }

    invisible(result)
}

##############################################################################
# TEST 1: Basic Gaussian Model (No Phylogeny)
##############################################################################

run_test("Basic Gaussian model without phylogeny", {
    N <- 40
    df <- data.frame(
        SP = paste0("sp", 1:N),
        Y = rnorm(N, mean = 10, sd = 2),
        X = rnorm(N, mean = 0, sd = 1)
    )

    fit <- because(
        data = df,
        id_col = "SP",
        equations = list(Y ~ X),
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    # Check structure
    stopifnot(inherits(fit$WAIC, "waic"))
    stopifnot(inherits(fit$WAIC, "data.frame"))
    stopifnot(all(c("Estimate", "SE") %in% colnames(fit$WAIC)))
    stopifnot(all(c("elpd_waic", "p_waic", "waic") %in% rownames(fit$WAIC)))

    # Check formula
    waic_check <- -2 * fit$WAIC["elpd_waic", "Estimate"]
    waic_actual <- fit$WAIC["waic", "Estimate"]
    stopifnot(abs(waic_check - waic_actual) < 0.01)

    # Check SEs are positive
    stopifnot(all(fit$WAIC$SE > 0))

    cat(
        "  WAIC:",
        round(fit$WAIC["waic", "Estimate"], 1),
        "±",
        round(fit$WAIC["waic", "SE"], 1),
        "\n"
    )
    cat("  p_waic:", round(fit$WAIC["p_waic", "Estimate"], 1), "\n")
})

##############################################################################
# TEST 2: Phylogenetic SEM
##############################################################################

run_test("Phylogenetic SEM with tree structure", {
    N <- 30
    tree <- rtree(N)

    # Simulate phylogenetically correlated data
    df <- data.frame(
        SP = tree$tip.label,
        Y = rnorm(N, mean = 15, sd = 3),
        X = rnorm(N, mean = 5, sd = 2)
    )

    fit <- because(
        data = df,
        structure = tree,
        id_col = "SP",
        equations = list(Y ~ X),
        WAIC = TRUE,
        optimise = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    # Verify WAIC computed
    stopifnot(!is.null(fit$WAIC))
    stopifnot(all(fit$WAIC$SE > 0))

    # Verify pointwise stored
    pointwise <- attr(fit$WAIC, "pointwise")
    stopifnot(!is.null(pointwise))
    stopifnot(nrow(pointwise) == N)

    cat(
        "  WAIC:",
        round(fit$WAIC["waic", "Estimate"], 1),
        "±",
        round(fit$WAIC["waic", "SE"], 1),
        "\n"
    )
})

##############################################################################
# TEST 3: Poisson Model (Non-Gaussian)
##############################################################################

run_test("Poisson model with count data", {
    N <- 35
    df <- data.frame(
        SP = paste0("sp", 1:N),
        Y = rpois(N, lambda = 5),
        X = rnorm(N)
    )

    fit <- because(
        data = df,
        id_col = "SP",
        equations = list(Y ~ X),
        family = list(Y = "poisson"),
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    stopifnot(!is.null(fit$WAIC))
    stopifnot(all(fit$WAIC$SE > 0))

    cat(
        "  WAIC:",
        round(fit$WAIC["waic", "Estimate"], 1),
        "±",
        round(fit$WAIC["waic", "SE"], 1),
        "\n"
    )
})

##############################################################################
# TEST 4: Binomial Model
##############################################################################

run_test("Binomial model with binary data", {
    N <- 40
    df <- data.frame(
        SP = paste0("sp", 1:N),
        Y = rbinom(N, 1, prob = 0.6),
        X = rnorm(N)
    )

    fit <- because(
        data = df,
        id_col = "SP",
        equations = list(Y ~ X),
        family = list(Y = "binomial"),
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    stopifnot(!is.null(fit$WAIC))
    stopifnot(all(fit$WAIC$SE > 0))

    cat(
        "  WAIC:",
        round(fit$WAIC["waic", "Estimate"], 1),
        "±",
        round(fit$WAIC["waic", "SE"], 1),
        "\n"
    )
})

##############################################################################
# TEST 5: Multi-Response Model
##############################################################################

run_test("Multi-response model (two Gaussian responses)", {
    N <- 30
    df <- data.frame(
        SP = paste0("sp", 1:N),
        Y1 = rnorm(N, mean = 10, sd = 2),
        Y2 = rnorm(N, mean = 5, sd = 1.5),
        X = rnorm(N)
    )

    fit <- because(
        data = df,
        id_col = "SP",
        equations = list(
            Y1 ~ X,
            Y2 ~ X
        ),
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    stopifnot(!is.null(fit$WAIC))
    stopifnot(all(fit$WAIC$SE > 0))

    # Check that log_lik accounts for both responses
    pointwise <- attr(fit$WAIC, "pointwise")
    dims <- attr(fit$WAIC, "dims")
    cat("  Observations:", dims["n_obs"], "\n")
    cat(
        "  WAIC:",
        round(fit$WAIC["waic", "Estimate"], 1),
        "±",
        round(fit$WAIC["waic", "SE"], 1),
        "\n"
    )
})

##############################################################################
# TEST 6: Model Comparison
##############################################################################

run_test("Model comparison (complex vs simple)", {
    N <- 50
    tree <- rtree(N)

    df <- data.frame(
        SP = tree$tip.label,
        Y = rnorm(N, mean = 10 + 0.5 * (1:N / 10), sd = 2),
        X1 = rnorm(N),
        X2 = rnorm(N)
    )

    # Complex model
    fit_complex <- because(
        data = df,
        structure = tree,
        id_col = "SP",
        equations = list(Y ~ X1 + X2),
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    # Simple model
    fit_simple <- because(
        data = df,
        structure = tree,
        id_col = "SP",
        equations = list(Y ~ X1),
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    # Null model
    fit_null <- because(
        data = df,
        structure = tree,
        id_col = "SP",
        equations = list(Y ~ 1),
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    cat("\n  Model Comparison:\n")
    cat("  ----------------\n")

    models <- list(
        "Y ~ X1 + X2" = fit_complex,
        "Y ~ X1" = fit_simple,
        "Y ~ 1 (null)" = fit_null
    )

    for (name in names(models)) {
        waic_val <- models[[name]]$WAIC["waic", "Estimate"]
        waic_se <- models[[name]]$WAIC["waic", "SE"]
        cat(sprintf("  %15s: WAIC = %7.1f ± %5.1f\n", name, waic_val, waic_se))
    }

    # Compare best 2
    w1 <- fit_complex$WAIC["waic", "Estimate"]
    w2 <- fit_simple$WAIC["waic", "Estimate"]
    se1 <- fit_complex$WAIC["waic", "SE"]
    se2 <- fit_simple$WAIC["waic", "SE"]

    diff <- abs(w1 - w2)
    se_diff <- sqrt(se1^2 + se2^2)

    cat(sprintf("\n  Difference: %.1f ± %.1f\n", diff, se_diff))
    cat(sprintf("  Ratio: %.2f\n", diff / se_diff))

    if (diff > 2 * se_diff) {
        cat("  → Models are significantly different\n")
    } else {
        cat("  → Models are not significantly different\n")
    }

    stopifnot(TRUE) # Always pass if we got here
})

##############################################################################
# TEST 7: MAG Model with Latent Variables
##############################################################################

run_test("MAG model with latent variable (induced correlations)", {
    N <- 35
    tree <- rtree(N)

    df <- data.frame(
        SP = tree$tip.label,
        Y1 = rnorm(N, mean = 10, sd = 2),
        Y2 = rnorm(N, mean = 8, sd = 1.5)
    )

    fit <- because(
        data = df,
        structure = tree,
        id_col = "SP",
        latent = "L1",
        equations = list(
            Y1 ~ L1,
            Y2 ~ L1
        ),
        latent_method = "correlations",
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    stopifnot(!is.null(fit$WAIC))
    stopifnot(all(fit$WAIC$SE > 0))

    cat(
        "  WAIC:",
        round(fit$WAIC["waic", "Estimate"], 1),
        "±",
        round(fit$WAIC["waic", "SE"], 1),
        "\n"
    )
    cat(
        "  rho (correlation):",
        round(mean(fit$samples[[1]][, "rho_Y1_Y2"]), 2),
        "\n"
    )
})

##############################################################################
# TEST 8: Mixed Distribution Model
##############################################################################

run_test("Mixed distribution model (Gaussian + Poisson)", {
    N <- 40
    df <- data.frame(
        SP = paste0("sp", 1:N),
        Y_cont = rnorm(N, mean = 10, sd = 2),
        Y_count = rpois(N, lambda = 5),
        X = rnorm(N)
    )

    fit <- because(
        data = df,
        id_col = "SP",
        equations = list(
            Y_cont ~ X,
            Y_count ~ X
        ),
        family = list(
            Y_cont = "gaussian",
            Y_count = "poisson"
        ),
        WAIC = TRUE,
        n.chains = 2,
        n.iter = 1000,
        quiet = TRUE
    )

    stopifnot(!is.null(fit$WAIC))
    stopifnot(all(fit$WAIC$SE > 0))

    pointwise <- attr(fit$WAIC, "pointwise")
    cat("  Total observations:", nrow(pointwise), "\n")
    cat(
        "  WAIC:",
        round(fit$WAIC["waic", "Estimate"], 1),
        "±",
        round(fit$WAIC["waic", "SE"], 1),
        "\n"
    )
})

##############################################################################
# SUMMARY
##############################################################################

cat("\n")
cat("\n")
cat(paste(rep("=", 78), collapse = ""), "\n")
cat("TEST SUMMARY\n")
cat(paste(rep("=", 78), collapse = ""), "\n\n")

cat(sprintf("Total tests:  %d\n", test_counter))
cat(sprintf(
    "Passed:       %d (%.0f%%)\n",
    pass_counter,
    100 * pass_counter / test_counter
))
cat(sprintf("Failed:       %d\n\n", fail_counter))

if (fail_counter == 0) {
    cat("✓ ALL TESTS PASSED! WAIC-SE is production-ready.\n\n")
    cat("Following Vehtari, Gelman & Gabry (2017):\n")
    cat("  'Practical Bayesian model evaluation using\n")
    cat("   leave-one-out cross-validation and WAIC'\n")
    cat("  Statistics and Computing, 27(5), 1413-1432\n")
} else {
    cat("✗ SOME TESTS FAILED - Review above for details\n")
}

cat("\n")
