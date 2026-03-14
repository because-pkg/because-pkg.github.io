# Test script for Ordinal Optimization in because
# Verifies correctness of the random effects formulation for ordinal models (Cumulative Logit)

# Source local files

library(rjags)
library(coda)
library(ape)
library(MASS)

set.seed(456)

# 1. Simulate Ordinal Data (e.g., IUCN Red List: LC=1, NT=2, VU=3, EN=4, CR=5)
N <- 100
tree <- rcoal(N)
VCV <- vcv(tree)
# Standardize VCV to max height 1
VCV <- VCV / max(node.depth.edgelength(tree))

# Predictor (e.g., log body mass)
X <- rnorm(N)

# True parameters
beta_true <- -0.8 # Larger body mass -> lower threat
lambda_true <- 0.6

# Phylogenetic + residual error
sigma_phylo <- sqrt(lambda_true)
sigma_resid <- sqrt(1 - lambda_true)
u <- mvrnorm(1, rep(0, N), sigma_phylo^2 * VCV)
epsilon <- rnorm(N, 0, sigma_resid)
eta <- beta_true * X + u + epsilon

# Cutpoints (true values)
cutpoints_true <- c(-1, 0, 1, 2)

# Calculate probabilities using cumulative logit
Q <- matrix(0, nrow = N, ncol = 4)
for (i in 1:N) {
    for (k in 1:4) {
        Q[i, k] <- 1 / (1 + exp(-(cutpoints_true[k] - eta[i])))
    }
}

# Category probabilities
P <- matrix(0, nrow = N, ncol = 5)
P[, 1] <- Q[, 1]
for (k in 2:4) {
    P[, k] <- Q[, k] - Q[, k - 1]
}
P[, 5] <- 1 - Q[, 4]

# Sample ordinal outcome
Y <- apply(P, 1, function(p) sample(1:5, 1, prob = p))

data <- list(Y = Y, X = X, K_Y = 5)
equations <- list(Y ~ X)

cat("Running Optimized Ordinal Model (Random Effects)...\n")
start_time <- Sys.time()
fit_opt <- because(
    data = data,
    tree = tree,
    equations = equations,
    family = c(Y = "ordinal"),
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
beta_est <- sum_opt["beta_Y_X", "Mean"]
lambda_est <- sum_opt["lambda_Y", "Mean"]

cat("\nComparison with True Values:\n")
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

# Check cutpoints
cat("\nEstimated Cutpoints:\n")
for (k in 1:4) {
    cutpoint_name <- paste0("cutpoint_Y[", k, "]")
    if (cutpoint_name %in% rownames(sum_opt)) {
        cat(
            "Cutpoint[",
            k,
            "]: True =",
            round(cutpoints_true[k], 3),
            "Estimated =",
            round(sum_opt[cutpoint_name, "Mean"], 3),
            "\n"
        )
    }
}

if (abs(beta_est - beta_true) < 0.3 && abs(lambda_est - lambda_true) < 0.3) {
    cat("\n✓ SUCCESS: Parameter recovery looks good!\n")
} else {
    cat("\n⚠ Parameter recovery could be improved (consider more iterations)\n")
}
