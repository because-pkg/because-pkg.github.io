# ==============================================================================
# Test Script: JAGS Optimization for Binomial Distribution
# ==============================================================================
# Verifies that the optimized random effects formulation works correctly
# for binomial (binary) response variables.
# ==============================================================================

# library(because)
library(ape)
library(coda)
library(rjags)

# Source local files

cat("\n=== Testing JAGS Optimization (Binomial Model) ===\n\n")

# 1. Simulate Data
set.seed(456)
N <- 50
tree <- rtree(N)
tree$edge.length <- tree$edge.length / max(branching.times(tree))
VCV <- vcv.phylo(tree)

# True parameters
beta <- 1.0 # Effect of X on logit(p)
lambda <- 0.6 # Phylogenetic signal
sigma_total <- 1.5 # Total variance on logit scale

# Calculate variance components
# lambda = phylo_var / (phylo_var + resid_var)
# total_var = phylo_var + resid_var
phylo_var <- lambda * sigma_total^2
resid_var <- (1 - lambda) * sigma_total^2

# Simulate X (predictor)
X <- rnorm(N)

# Simulate phylogenetic + residual error on logit scale
Sigma_phylo <- phylo_var * VCV
err_phylo <- MASS::mvrnorm(1, rep(0, N), Sigma_phylo)
err_resid <- rnorm(N, 0, sqrt(resid_var))
err_total <- err_phylo + err_resid

# Calculate probability and simulate binary outcomes
logit_p <- beta * X + err_total
p <- plogis(logit_p)
Y <- rbinom(N, 1, p)

cat(sprintf(
    "Simulated data: %d binary observations (%.1f%% ones)\n",
    N,
    100 * mean(Y)
))

data <- list(Y = Y, X = X)
equations <- list(Y ~ X)

# 2. Run Optimized Model
cat("\nRunning optimized model (optimize = TRUE)...\n")
time_opt <- system.time({
    fit_opt <- because(
        data = data,
        tree = tree,
        equations = equations,
        family = c(Y = "binomial"),
        n.iter = 3000,
        n.burnin = 500,
        n.thin = 3,
        optimise = TRUE,
        quiet = TRUE
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
cat("\nChecking convergence (R-hat):\n")
gelman <- gelman.diag(fit_opt$samples)
print(gelman)

# Focus on key parameters
key_params <- c("beta_Y_X", "lambdaY")
if (all(key_params %in% rownames(gelman$psrf))) {
    if (all(gelman$psrf[key_params, "Point est."] < 1.1)) {
        cat("✓ Convergence successful for key parameters (R-hat < 1.1)\n")
    } else {
        cat("⚠ Convergence warning for key parameters\n")
    }
}

# Check parameter recovery
if ("beta_Y_X" %in% rownames(sum_opt$statistics)) {
    beta_est <- sum_opt$statistics["beta_Y_X", "Mean"]
    cat(sprintf("\nTrue beta: %.2f, Estimated: %.2f\n", beta, beta_est))
}

if ("lambdaY" %in% rownames(sum_opt$statistics)) {
    lambda_est <- sum_opt$statistics["lambdaY", "Mean"]
    cat(sprintf("True lambda: %.2f, Estimated: %.2f\n", lambda, lambda_est))

    if (abs(lambda_est - lambda) < 0.3) {
        cat("✓ Lambda recovery reasonable (within 0.3)\n")
    } else {
        cat("⚠ Lambda estimate differs from true value\n")
    }
}

# 4. Run Unoptimized Model (for comparison)
cat("\nRunning unoptimized model (optimize = FALSE)...\n")
time_unopt <- system.time({
    fit_unopt <- because(
        data = data,
        tree = tree,
        equations = equations,
        family = c(Y = "binomial"),
        n.iter = 3000,
        n.burnin = 500,
        n.thin = 3,
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
if (
    "beta_Y_X" %in%
        rownames(sum_unopt$statistics) &&
        "beta_Y_X" %in% rownames(sum_opt$statistics)
) {
    beta_unopt <- sum_unopt$statistics["beta_Y_X", "Mean"]
    beta_opt <- sum_opt$statistics["beta_Y_X", "Mean"]
    cat(sprintf(
        "\nOptimized beta: %.4f, Unoptimized: %.4f\n",
        beta_opt,
        beta_unopt
    ))
}

if (
    "lambdaY" %in%
        rownames(sum_unopt$statistics) &&
        "lambdaY" %in% rownames(sum_opt$statistics)
) {
    lambda_unopt <- sum_unopt$statistics["lambdaY", "Mean"]
    lambda_opt <- sum_opt$statistics["lambdaY", "Mean"]
    cat(sprintf(
        "Optimized lambda: %.4f, Unoptimized: %.4f\n",
        lambda_opt,
        lambda_unopt
    ))

    if (abs(lambda_opt - lambda_unopt) < 0.1) {
        cat("✓ Estimates match between methods (diff < 0.1)\n")
    } else {
        cat("⚠ Estimates differ between methods\n")
    }
}

cat("\n=== Binomial Optimization Test Complete ===\n")
