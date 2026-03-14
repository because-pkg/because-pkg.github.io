# ==============================================================================
# Test Script: JAGS Optimization (Simple Model)
# ==============================================================================
# Verifies that the optimized random effects formulation works correctly
# for a simple regression model.
# ==============================================================================

# library(because)
library(ape)
library(coda)
library(rjags)

# Source local files

cat("\n=== Testing JAGS Optimization (Simple Model) ===\n\n")

# 1. Simulate Data
set.seed(123)
N <- 50
tree <- rtree(N)
tree$edge.length <- tree$edge.length / max(branching.times(tree))
VCV <- vcv.phylo(tree)

# True parameters
beta <- 0.5
lambda <- 0.7
sigma <- 1.0
tau <- 1 / sigma^2

# Simulate Y ~ N(beta*X, sigma^2 * (lambda*V + (1-lambda)*I))
X <- rnorm(N)
mu <- beta * X
Sigma <- sigma^2 * (lambda * VCV + (1 - lambda) * diag(N))
Y <- MASS::mvrnorm(1, mu, Sigma)

data <- list(Y = Y, X = X)
equations <- list(Y ~ X)

# Debug: Print model code
cat("Generating model code...\n")
model_out <- because_model(
    equations = equations,
    optimise = TRUE
)
cat("---------------------------------------------------\n")
cat(model_out$model)
cat("\n---------------------------------------------------\n")

# 2. Run Optimized Model
cat("Running optimized model (optimize = TRUE)...\n")
time_opt <- system.time({
    fit_opt <- because(
        data = data,
        tree = tree,
        equations = equations,
        n.iter = 5000,
        n.burnin = 1000,
        n.thin = 5,
        optimise = TRUE,
        quiet = FALSE # Enable output to see errors
    )
})

cat(sprintf(
    "✓ Optimized model finished in %.2f seconds\n",
    time_opt["elapsed"]
))

# 3. Verify Output
cat("\nChecking parameter estimates:\n")
sum_opt <- fit_opt$summary
print(sum_opt$statistics[, c("Mean", "SD")])

# Check convergence
# Check convergence
cat("\nChecking convergence (R-hat):\n")
gelman <- gelman.diag(fit_opt$samples)
print(gelman)

if (all(gelman$psrf[c("beta_Y_X", "lambdaY"), "Point est."] < 1.1)) {
    cat("✓ Convergence successful for key parameters (R-hat < 1.1)\n")
} else {
    cat("⚠ Convergence warning for key parameters\n")
}

# Check parameter recovery
beta_est <- sum_opt$statistics["beta_Y_X", "Mean"]
lambda_est <- sum_opt$statistics["lambdaY", "Mean"]

cat(sprintf("\nTrue beta: %.2f, Estimated: %.2f\n", beta, beta_est))
cat(sprintf("True lambda: %.2f, Estimated: %.2f\n", lambda, lambda_est))

if (abs(beta_est - beta) < 0.2 && abs(lambda_est - lambda) < 0.2) {
    cat("✓ Parameter recovery successful\n")
} else {
    cat("⚠ Parameter recovery warning\n")
}

# 4. Run Unoptimized Model (for comparison)
cat("\nRunning unoptimized model (optimize = FALSE)...\n")
time_unopt <- system.time({
    fit_unopt <- because(
        data = data,
        tree = tree,
        equations = equations,
        n.iter = 5000,
        n.burnin = 1000,
        n.thin = 5,
        optimise = FALSE,
        quiet = TRUE
    )
})

cat(sprintf(
    "✓ Unoptimized model finished in %.2f seconds\n",
    time_unopt["elapsed"]
))
cat(sprintf("Speedup: %.2fx\n", time_unopt["elapsed"] / time_opt["elapsed"]))

# Compare estimates
sum_unopt <- fit_unopt$summary
beta_unopt <- sum_unopt$statistics["beta_Y_X", "Mean"]
lambda_unopt <- sum_unopt$statistics["lambdaY", "Mean"]

cat(sprintf(
    "\nOptimized beta: %.4f, Unoptimized: %.4f\n",
    beta_est,
    beta_unopt
))
cat(sprintf(
    "Optimized lambda: %.4f, Unoptimized: %.4f\n",
    lambda_est,
    lambda_unopt
))

if (abs(beta_est - beta_unopt) < 0.05) {
    cat("✓ Estimates match between methods\n")
} else {
    cat("⚠ Estimates differ between methods\n")
}
