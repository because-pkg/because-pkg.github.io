# Test script for Multinomial Optimization in because
# Verifies correctness and speedup of the random effects formulation for multinomial models

# Source local files

library(rjags)
library(coda)
library(ape)
library(MASS)

set.seed(123)

# 1. Simulate Data
N <- 100
tree <- rcoal(N)
VCV <- vcv(tree)
# Standardize VCV to max height 1 (matches because logic)
VCV <- VCV / max(node.depth.edgelength(tree))

# Predictor
X <- rnorm(N)
beta_true <- c(0, 0.5, -0.5) # k=1 (ref), k=2, k=3
alpha_true <- c(0, -1, 1)

# Phylogenetic signal (lambda = 0.6)
lambda_true <- 0.6
sigma_phylo <- sqrt(lambda_true)
sigma_resid <- sqrt(1 - lambda_true)

# Latent variables for k=2, 3 (k=1 is reference 0)
L <- matrix(0, nrow = N, ncol = 3)
# k=1 is 0
L[, 1] <- 0

# k=2
u2 <- mvrnorm(1, rep(0, N), sigma_phylo^2 * VCV)
e2 <- rnorm(N, 0, sigma_resid)
L[, 2] <- alpha_true[2] + beta_true[2] * X + u2 + e2

# k=3
u3 <- mvrnorm(1, rep(0, N), sigma_phylo^2 * VCV)
e3 <- rnorm(N, 0, sigma_resid)
L[, 3] <- alpha_true[3] + beta_true[3] * X + u3 + e3

# Probabilities (Softmax)
P <- exp(L) / rowSums(exp(L))

# Sample Y (Categorical 1, 2, 3)
Y <- apply(P, 1, function(p) sample(1:3, 1, prob = p))

data <- list(Y = Y, X = X)
equations <- list(Y ~ X)

cat("Running Unoptimized Multinomial Model (Marginal)...\n")
start_time <- Sys.time()
fit_marg <- because(
    data = data,
    tree = tree,
    equations = equations,
    family = c(Y = "multinomial"),
    optimise = FALSE,
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 2 # Reduced for speed
)
end_time <- Sys.time()
time_marg <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Unoptimized Time:", time_marg, "s\n")

cat("\nRunning Optimized Multinomial Model (Random Effects)...\n")
start_time <- Sys.time()
fit_opt <- because(
    data = data,
    tree = tree,
    equations = equations,
    family = c(Y = "multinomial"),
    optimise = TRUE,
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 2
)
end_time <- Sys.time()
time_opt <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Optimized Time:", time_opt, "s\n")

speedup <- time_marg / time_opt
cat("\nSpeedup:", round(speedup, 2), "x\n")

# 3. Verify Results
sum_opt <- fit_opt$summary$statistics
sum_marg <- fit_marg$summary$statistics

print(head(sum_opt))

# Check convergence (Rhat)
max_rhat <- max(sum_opt[, "Rhat"], na.rm = TRUE)
cat("\nMax R-hat (Optimized):", max_rhat, "\n")

if (max_rhat > 1.1) {
    warning("Optimized model has convergence issues (Rhat > 1.1)")
}

# Compare estimates for beta_Y_X[2] and beta_Y_X[3]
# Note: Parameter names might be indexed like beta_Y_X[2]
beta2_opt <- sum_opt["beta_Y_X[2]", "Mean"]
beta2_marg <- sum_marg["beta_Y_X[2]", "Mean"]
beta3_opt <- sum_opt["beta_Y_X[3]", "Mean"]
beta3_marg <- sum_marg["beta_Y_X[3]", "Mean"]

cat("\nComparison of Estimates:\n")
cat(
    "Beta[2]: Optimized =",
    round(beta2_opt, 3),
    "Marginal =",
    round(beta2_marg, 3),
    "\n"
)
cat(
    "Beta[3]: Optimized =",
    round(beta3_opt, 3),
    "Marginal =",
    round(beta3_marg, 3),
    "\n"
)

# Check Lambda estimates
lambda2_opt <- sum_opt["lambda_Y[2]", "Mean"]
lambda2_marg <- sum_marg["lambda_Y[2]", "Mean"]
cat(
    "Lambda[2]: Optimized =",
    round(lambda2_opt, 3),
    "Marginal =",
    round(lambda2_marg, 3),
    "\n"
)

if (abs(beta2_opt - beta2_marg) > 0.5) {
    warning("Significant difference in Beta[2] estimates!")
} else {
    cat("\nSUCCESS: Estimates match reasonably well.\n")
}
