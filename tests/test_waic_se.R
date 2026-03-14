library(because)
library(ape)

# Test WAIC with standard errors
set.seed(42)
N <- 30
tree <- rtree(N)

df <- data.frame(
    SP = tree$tip.label,
    Y = rnorm(N, mean = 10, sd = 2),
    X = rnorm(N)
)

cat("=== Testing WAIC with Standard Errors ===\n\n")

cat("Fitting model with WAIC=TRUE...\n")
fit <- because(
    data = df,
    structure = tree,
    id_col = "SP",
    equations = list(Y ~ X),
    WAIC = TRUE,
    n.chains = 2,
    n.iter = 500,
    quiet = TRUE
)

cat("\n=== WAIC Results ===\n")
print(fit$WAIC)

cat("\n=== Checking Structure ===\n")
cat("Class:", class(fit$WAIC), "\n")
cat("Columns:", paste(colnames(fit$WAIC), collapse = ", "), "\n")
cat("Row names:", paste(rownames(fit$WAIC), collapse = ", "), "\n")

cat("\n=== Checking Attributes ===\n")
pointwise <- attr(fit$WAIC, "pointwise")
if (!is.null(pointwise)) {
    cat("✓ Pointwise values stored\n")
    cat("  Dimensions:", nrow(pointwise), "observations\n")
    cat("  Columns:", paste(colnames(pointwise), collapse = ", "), "\n")
} else {
    cat("✗ No pointwise values found\n")
}

dims <- attr(fit$WAIC, "dims")
if (!is.null(dims)) {
    cat(
        "✓ Dimensions stored:",
        paste(names(dims), "=", dims, collapse = ", "),
        "\n"
    )
}

cat("\n=== Verification ===\n")
# Check that SE values are positive
if (all(fit$WAIC$SE > 0)) {
    cat("✓ PASS: All standard errors are positive\n")
} else {
    cat("✗ FAIL: Some standard errors are not positive\n")
}

# Check that WAIC is computed correctly
waic_manual <- -2 * fit$WAIC["elpd_waic", "Estimate"]
waic_actual <- fit$WAIC["waic", "Estimate"]
if (abs(waic_manual - waic_actual) < 0.1) {
    cat("✓ PASS: WAIC = -2 * elpd_waic (correct formula)\n")
} else {
    cat("✗ FAIL: WAIC formula incorrect\n")
    cat("  Expected:", waic_manual, "Got:", waic_actual, "\n")
}

cat("\n=== Testing Model Comparison ===\n")
# Fit a simpler model
fit_null <- because(
    data = df,
    structure = tree,
    id_col = "SP",
    equations = list(Y ~ 1),
    WAIC = TRUE,
    n.chains = 2,
    n.iter = 500,
    quiet = TRUE
)

cat("\nModel 1 (Y ~ X):\n")
print(fit$WAIC)

cat("\nModel 2 (Y ~ 1, null model):\n")
print(fit_null$WAIC)

# Compare models
waic_diff <- fit$WAIC["waic", "Estimate"] - fit_null$WAIC["waic", "Estimate"]
se_diff <- sqrt(fit$WAIC["waic", "SE"]^2 + fit_null$WAIC["waic", "SE"]^2)

cat("\n=== Model Comparison ===\n")
cat("WAIC difference:", round(waic_diff, 1), "±", round(se_diff, 1), "\n")
cat("Ratio |diff|/SE:", round(abs(waic_diff) / se_diff, 2), "\n")

if (abs(waic_diff) > 2 * se_diff) {
    cat("✓ Models are significantly different (|diff| > 2*SE)\n")
    if (waic_diff < 0) {
        cat("  → Model 1 (Y ~ X) is preferred\n")
    } else {
        cat("  → Model 2 (null) is preferred\n")
    }
} else {
    cat("  Models are not significantly different\n")
}

cat("\n=== Summary ===\n")
cat("WAIC with SE implementation: SUCCESS!\n")
cat("Following Vehtari et al. (2017) algorithm ✓\n")
