# Test script for Latent Variable implementation in because

library(because)
library(ape)

# 1. Simulate data with a latent common cause
# L -> X
# L -> Y
# This induces correlation between X and Y residuals

set.seed(123)
N <- 50
tree <- rtree(N)
tree$edge.length <- tree$edge.length / max(branching.times(tree))

# Latent variable L (phylogenetically structured)
L <- rTraitCont(tree, model = "BM", sigma = 1)

# X and Y depend on L plus independent error
X <- 0.5 * L + rnorm(N, 0, 0.5)
Y <- -0.5 * L + rnorm(N, 0, 0.5)

data <- list(
    X = X,
    Y = Y,
    N = N
)

# 2. Define model with latent variable
# We specify L in the equations but mark it as latent
equations <- list(
    X ~ L,
    Y ~ L
)

# 3. Run because with latent argument
# We expect:
# - d-sep to identify induced correlation between X and Y
# - JAGS model to include correlated residuals for X and Y

cat("Running because with latent variable...\n")
fit <- because(
    data = data,
    tree = tree,
    equations = equations,
    latent = "L",
    dsep = TRUE,
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 2,
    quiet = FALSE
)

# 4. Inspect results

# Check if induced correlations were found
if (!is.null(fit$induced_correlations)) {
    cat("\nSUCCESS: Induced correlations found:\n")
    print(fit$induced_correlations)
} else {
    cat("\nFAILURE: No induced correlations found.\n")
}

# Check JAGS model code for correlated residuals
# Check if model contains correlated residuals
model_code <- fit$model_code
if (is.null(model_code)) {
    cat("fit$model_code is NULL!\n")
} else if (!is.character(model_code)) {
    cat("fit$model_code is NOT character! Type:", typeof(model_code), "\n")
}

if (grepl("rho_X_Y", model_code) || grepl("rho_Y_X", model_code)) {
    cat(
        "\nSUCCESS: JAGS model contains rho parameter for correlated residuals.\n"
    )
} else {
    cat("\nFAILURE: JAGS model missing rho parameter.\n")
    cat("Model code snippet:\n")
    cat(substr(model_code, 1, 1000)) # Print first 1000 chars
}

# Check posterior for rho
cat("\nMonitored parameters:\n")
print(fit$monitor)
cat("\nSummary statistics rownames:\n")
print(rownames(fit$summary$statistics))
cat("\nSummary structure:\n")
str(fit$summary)
cat("\nSamples class:\n")
print(class(fit$samples))
cat("\nSamples length:\n")
print(length(fit$samples))
if (length(fit$samples) > 0) {
    cat("\nChain 1 dimensions:\n")
    print(dim(fit$samples[[1]]))
}

cat("\nTrying summary locally:\n")
library(coda)
local_summary <- summary(fit$samples)
print(str(local_summary))

if (
    "rho_X_Y" %in%
        rownames(fit$summary$statistics) ||
        "rho_Y_X" %in% rownames(fit$summary$statistics)
) {
    cat("\nSUCCESS: rho parameter estimated in summary.\n")
    print(fit$summary$statistics[
        grep("rho", rownames(fit$summary$statistics)),
    ])
} else {
    cat("\nFAILURE: rho parameter not in summary.\n")
}

cat("\nTest complete.\n")
