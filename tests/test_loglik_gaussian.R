library(because)

# Create simple test data
set.seed(123)
N <- 20
df <- data.frame(
    SP = paste0("sp", 1:N),
    Y = rnorm(N, mean = 10, sd = 2),
    X = rnorm(N)
)

# Run model and check if log_lik is in the JAGS model
cat("Running simple Gaussian model...\n")
fit <- because(
    data = df,
    id_col = "SP",
    equations = list(Y ~ X),
    n.chains = 1,
    n.iter = 100,
    quiet = FALSE
)

# Check if model code contains log_lik
cat("\n=== Checking JAGS Model for log_lik ===\n")
model_str <- fit$model_code

if (grepl("log_lik", model_str)) {
    cat("✓ PASS: log_lik found in model\n")

    # Extract log_lik lines
    model_lines <- strsplit(model_str, "\n")[[1]]
    log_lik_lines <- grep("log_lik", model_lines, value = TRUE)
    cat("\nlog_lik statements:\n")
    cat(paste(log_lik_lines, collapse = "\n"), "\n")
} else {
    cat("✗ FAIL: log_lik NOT found in model\n")
}

# Check if log_lik would be in monitored parameters if we had enabled it
cat("\n=== Summary ===\n")
cat("Model compiled successfully!\n")
cat("Next step: Add log_lik to monitored parameters in because()\n")
