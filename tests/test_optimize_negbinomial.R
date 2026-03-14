# Test script for Negative Binomial Optimization in because
# Verifies correctness of the random effects formulation for NB models

# Source local files

library(rjags)
library(coda)
library(ape)
library(MASS)

set.seed(321)

# 1. Simulate Negative Binomial Data (heavily overdispersed)
N <- 100
tree <- rcoal(N)
VCV <- vcv(tree)
# Standardize VCV to max height 1
VCV <- VCV / max(node.depth.edgelength(tree))

# Predictor (e.g., habitat quality)
X <- rnorm(N)

# True parameters
beta_true <- 0.6 # Positive effect
alpha_true <- 3 # Baseline log(count)
lambda_true <- 0.4
r_true <- 2 # Size parameter (smaller = more overdispersion)

# Phylogenetic + residual error
sigma_phylo <- sqrt(lambda_true)
sigma_resid <- sqrt(1 - lambda_true)
u <- mvrnorm(1, rep(0, N), sigma_phylo^2 * VCV)
epsilon <- rnorm(N, 0, sigma_resid)

# Log(mean)
log_mu <- alpha_true + beta_true * X + u + epsilon
mu <- exp(log_mu)

# Sample NB counts
# In R: rnbinom(n, size, mu) - size is dispersion, lower = more overdispersed
Y <- rnbinom(N, size = r_true, mu = mu)

data <- list(Y = Y, X = X)
equations <- list(Y ~ X)

cat("Running Optimized Negative Binomial Model...\n")
cat("True overdispersion (r):", r_true, "\n")
cat("(Smaller r = more overdispersion)\n\n")

start_time <- Sys.time()
fit_opt <- because(
    data = data,
    tree = tree,
    equations = equations,
    family = c(Y = "negbinomial"),
    optimise = TRUE,
    n.iter = 2000,
    n.burnin = 1000,
    n.chains = 2
)
end_time <- Sys.time()
time_opt <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Optimized Time:", time_opt, "s\n")

# 3. Verify Results
sum_opt <- fit_opt$summary$statistics

cat("\nParameter Estimates (First 12):\n")
print(head(sum_opt, 12))

# Check convergence (Rhat)
if ("Rhat" %in% colnames(sum_opt)) {
    max_rhat <- max(sum_opt[, "Rhat"], na.rm = TRUE)
    cat("\nMax R-hat (Optimized):", max_rhat, "\n")

    if (max_rhat > 1.1) {
        warning("Optimized model has convergence issues (Rhat > 1.1)")
    } else {
        cat("✓ Convergence successful (R-hat < 1.1)\n")
    }
}

# Check parameter recovery
alpha_est <- sum_opt["alphaY", "Mean"]
beta_est <- sum_opt["beta_Y_X", "Mean"]
lambda_est <- sum_opt["lambda_Y", "Mean"]

cat("\nComparison with True Values:\n")
cat(
    "Alpha: True =",
    round(alpha_true, 3),
    "Estimated =",
    round(alpha_est, 3),
    "\n"
)
cat(
    "Beta: True =",
    round(beta_true, 3),
    "Estimated =",
    round(beta_est, 3),
    "\n"
)
cat(
    "Lambda: True =",
    round(lambda_true, 3),
    "Estimated =",
    round(lambda_est, 3),
    "\n"
)

if ("r_Y" %in% rownames(sum_opt)) {
    r_est <- sum_opt["r_Y", "Mean"]
    cat(
        "Size (r): True =",
        round(r_true, 3),
        "Estimated =",
        round(r_est, 3),
        "\n"
    )
    cat("\nOverdispersion interpretation:\n")
    cat("Estimated r =", round(r_est, 2), "(smaller = more overdispersion)\n")
} else {
    cat("Note: r_Y not in summary (add 'r' to monitor list)\n")
}

cat(
    "Variance/Mean ratio (theoretical) =",
    round(1 + mean(mu) / r_true, 2),
    "\n"
)

if (abs(beta_est - beta_true) < 0.4 && abs(lambda_est - lambda_true) < 0.3) {
    cat("\n✓ SUCCESS: Parameter recovery looks good!\n")
} else {
    cat("\n⚠ Parameter recovery could be improved (consider more iterations)\n")
}
