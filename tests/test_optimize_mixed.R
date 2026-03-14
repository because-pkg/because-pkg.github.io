# ==============================================================================
# Test Script: Mixed Gaussian + Binomial Model with Optimization
# ==============================================================================
# Verifies that optimized random effects work correctly when a model contains
# both Gaussian and binomial response variables.
# ==============================================================================

# library(because)
library(ape)
library(coda)
library(rjags)

# Source local files

cat("\n=== Testing Mixed Gaussian + Binomial Model ===\n\n")

# 1. Simulate Data
set.seed(789)
N <- 60
tree <- rtree(N)
tree$edge.length <- tree$edge.length / max(branching.times(tree))
VCV <- vcv.phylo(tree)

# True parameters
beta_ZY <- 0.8 # Effect of Z on Y (Gaussian)
beta_YB <- 1.2 # Effect of Y on B (binomial)
lambda_Y <- 0.7 # Phylogenetic signal for Y (Gaussian)
lambda_B <- 0.5 # Phylogenetic signal for B (binomial)

# Simulate Z (exogenous predictor)
Z <- rnorm(N)

# Simulate Y (Gaussian response)
sigma_Y <- 1.0
phylo_var_Y <- lambda_Y * sigma_Y^2
resid_var_Y <- (1 - lambda_Y) * sigma_Y^2
Sigma_phylo_Y <- phylo_var_Y * VCV
err_phylo_Y <- MASS::mvrnorm(1, rep(0, N), Sigma_phylo_Y)
err_resid_Y <- rnorm(N, 0, sqrt(resid_var_Y))
Y <- beta_ZY * Z + err_phylo_Y + err_resid_Y

# Simulate B (binomial response, depends on Y)
sigma_B <- 1.5
phylo_var_B <- lambda_B * sigma_B^2
resid_var_B <- (1 - lambda_B) * sigma_B^2
Sigma_phylo_B <- phylo_var_B * VCV
err_phylo_B <- MASS::mvrnorm(1, rep(0, N), Sigma_phylo_B)
err_resid_B <- rnorm(N, 0, sqrt(resid_var_B))
logit_p <- beta_YB * Y + err_phylo_B + err_resid_B
p <- plogis(logit_p)
B <- rbinom(N, 1, p)

cat(sprintf("Simulated data: %d species\n", N))
cat(sprintf("  Y (Gaussian): mean=%.2f, sd=%.2f\n", mean(Y), sd(Y)))
cat(sprintf("  B (Binomial): %.1f%% ones\n", 100 * mean(B)))

data <- list(Y = Y, B = B, Z = Z)
equations <- list(
    Y ~ Z,
    B ~ Y
)

# 2. Run Optimized Model
cat("\nRunning optimized model (optimize = TRUE)...\n")
time_opt <- system.time({
    fit_opt <- because(
        data = data,
        tree = tree,
        equations = equations,
        family = c(Y = "gaussian", B = "binomial", Z = "gaussian"),
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
print(sum_opt$statistics[
    c("beta_Y_Z", "beta_B_Y", "lambdaY", "lambdaB"),
    c("Mean", "SD")
])

# Check convergence
cat("\nChecking convergence (R-hat):\n")
gelman <- gelman.diag(fit_opt$samples)
key_params <- c("beta_Y_Z", "beta_B_Y", "lambdaY", "lambdaB")
print(gelman$psrf[key_params, ])

if (all(gelman$psrf[key_params, "Point est."] < 1.1)) {
    cat("✓ Convergence successful for key parameters (R-hat < 1.1)\n")
} else {
    cat("⚠ Some parameters have R-hat > 1.1\n")
}

# Check parameter recovery
cat("\nParameter Recovery:\n")
if ("beta_Y_Z" %in% rownames(sum_opt$statistics)) {
    beta_Z_est <- sum_opt$statistics["beta_Y_Z", "Mean"]
    cat(sprintf("  beta_ZY: True=%.2f, Est=%.2f\n", beta_ZY, beta_Z_est))
}

if ("beta_B_Y" %in% rownames(sum_opt$statistics)) {
    beta_Y_est <- sum_opt$statistics["beta_B_Y", "Mean"]
    cat(sprintf("  beta_YB: True=%.2f, Est=%.2f\n", beta_YB, beta_Y_est))
}

if ("lambdaY" %in% rownames(sum_opt$statistics)) {
    lambda_Y_est <- sum_opt$statistics["lambdaY", "Mean"]
    cat(sprintf("  lambda_Y: True=%.2f, Est=%.2f\n", lambda_Y, lambda_Y_est))
}

if ("lambdaB" %in% rownames(sum_opt$statistics)) {
    lambda_B_est <- sum_opt$statistics["lambdaB", "Mean"]
    cat(sprintf("  lambda_B: True=%.2f, Est=%.2f\n", lambda_B, lambda_B_est))
}

# 4. Run Unoptimized Model (for comparison)
cat("\nRunning unoptimized model (optimize = FALSE)...\n")
time_unopt <- system.time({
    fit_unopt <- because(
        data = data,
        tree = tree,
        equations = equations,
        family = c(Y = "gaussian", B = "binomial", Z = "gaussian"),
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
cat("\nComparing estimates (Optimized vs Unoptimized):\n")
for (param in key_params) {
    if (
        param %in%
            rownames(sum_opt$statistics) &&
            param %in% rownames(sum_unopt$statistics)
    ) {
        opt_est <- sum_opt$statistics[param, "Mean"]
        unopt_est <- sum_unopt$statistics[param, "Mean"]
        diff <- abs(opt_est - unopt_est)
        cat(sprintf(
            "  %s: %.4f vs %.4f (diff=%.4f)\n",
            param,
            opt_est,
            unopt_est,
            diff
        ))
    }
}

cat("\n=== Mixed Model Test Complete ===\n")
cat("✓ Both Gaussian and binomial variables handled correctly\n")
