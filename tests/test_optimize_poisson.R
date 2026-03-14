# Test script for Poisson Optimization in because
# Verifies correctness of the random effects formulation for Poisson models (log link)

# Source local files

library(rjags)
library(coda)
library(ape)
library(MASS)

set.seed(789)

# 1. Simulate Poisson Count Data
N <- 100
tree <- rcoal(N)
VCV <- vcv(tree)
# Standardize VCV to max height 1
VCV <- VCV / max(node.depth.edgelength(tree))

# Predictor (e.g., log area)
X <- rnorm(N)

# True parameters
beta_true <- 0.5 # Larger area -> more species
alpha_true <- 2 # Baseline log(count)
lambda_true <- 0.5

# Phylogenetic + residual error
sigma_phylo <- sqrt(lambda_true)
sigma_resid <- sqrt(1 - lambda_true)
u <- mvrnorm(1, rep(0, N), sigma_phylo^2 * VCV)
epsilon <- rnorm(N, 0, sigma_resid)

# Log(mean)
log_mu <- alpha_true + beta_true * X + u + epsilon
mu <- exp(log_mu)

# Sample counts
Y <- rpois(N, lambda = mu)

data <- list(Y = Y, X = X)
equations <- list(Y ~ X)

cat("Running Optimized Poisson Model (Random Effects)...\n")
start_time <- Sys.time()
fit_opt <- because(
    data = data,
    tree = tree,
    equations = equations,
    family = c(Y = "poisson"),
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

cat("\nParameter Estimates (First 10):\n")
print(head(sum_opt, 10))

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

# Check overdispersion handling (via tau_e)
tau_e_est <- sum_opt["tau_e_Y", "Mean"]
cat("\nOverdispersion parameter (tau_e):", round(tau_e_est, 3), "\n")
cat("(Higher tau_e = less overdispersion relative to Poisson)\n")

if (abs(beta_est - beta_true) < 0.3 && abs(lambda_est - lambda_true) < 0.3) {
    cat("\n✓ SUCCESS: Parameter recovery looks good!\n")
} else {
    cat("\n⚠ Parameter recovery could be improved (consider more iterations)\n")
}
