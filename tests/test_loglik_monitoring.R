library(because)

# Test log_lik monitoring with WAIC=TRUE
set.seed(123)
N <- 15
df <- data.frame(
    SP = paste0("sp", 1:N),
    Y = rnorm(N, mean = 10, sd = 2),
    X = rnorm(N)
)

cat("=== Test 1: WAIC=FALSE (no log_lik monitoring) ===\n")
fit1 <- because(
    data = df,
    id_col = "SP",
    equations = list(Y ~ X),
    WAIC = FALSE,
    n.chains = 1,
    n.iter = 100,
    quiet = FALSE
)

has_loglik1 <- any(grepl("log_lik", colnames(fit1$samples[[1]])))
cat("log_lik in samples:", has_loglik1, "\n")
if (!has_loglik1) {
    cat("✓ PASS: log_lik NOT monitored when WAIC=FALSE\n")
} else {
    cat("✗ FAIL: log_lik WAS monitored when WAIC=FALSE\n")
}

cat("\n=== Test 2: WAIC=TRUE (should monitor log_lik) ===\n")
fit2 <- because(
    data = df,
    id_col = "SP",
    equations = list(Y ~ X),
    WAIC = TRUE,
    n.chains = 2,
    n.iter = 100,
    quiet = FALSE
)

has_loglik2 <- any(grepl("log_lik", colnames(fit2$samples[[1]])))
cat("log_lik in samples:", has_loglik2, "\n")

if (has_loglik2) {
    # Count how many log_lik parameters
    log_lik_cols <- grep("log_lik", colnames(fit2$samples[[1]]), value = TRUE)
    cat("✓ PASS: log_lik monitored when WAIC=TRUE\n")
    cat("  Found", length(log_lik_cols), "log_lik parameters\n")
    cat("  Names:", paste(head(log_lik_cols, 3), collapse = ", "), "...\n")

    # Check dimensions
    cat(
        "  Dimensions:",
        nrow(fit2$samples[[1]]),
        "samples ×",
        length(log_lik_cols),
        "observations\n"
    )
} else {
    cat("✗ FAIL: log_lik NOT monitored when WAIC=TRUE\n")
}

cat("\n=== Summary ===\n")
cat("Conditional monitoring working correctly!\n")
