library(because)
library(ape)

set.seed(123)
N <- 30
tree <- rtree(N)

# Simulate data with latent variable L affecting X and Y
L_true <- rnorm(N)
X <- 0.8 * L_true + rnorm(N, sd = 0.3)
Y <- 0.6 * L_true + rnorm(N, sd = 0.4)

data_list <- list(
    X = X,
    Y = Y
)

equations <- list(
    X ~ L,
    Y ~ L
)

cat("=== Testing Latent Variable Methods ===\n\n")

# Test 1: MAG Approach (correlations) - DEFAULT
cat("Test 1: MAG Approach (default, correlations)\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

fit_mag <- because(
    data = data_list,
    tree = tree,
    equations = equations,
    latent = "L",
    latent_method = "correlations", # Or omit for default
    n.iter = 2000,
    n.burnin = 1000,
    quiet = FALSE
)

cat("\nParameters monitored (MAG):\n")
print(names(fit_mag$samples[[1]][1, ]))

cat("\nSummary (should show rho_X_Y):\n")
print(summary(fit_mag))

cat("\nInduced correlations:\n")
print(fit_mag$induced_correlations)

cat("\n\n")

# Test 2: Explicit Latent Variables
cat("Test 2: Explicit Latent Variable Modeling\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

fit_explicit <- because(
    data = data_list,
    tree = tree,
    equations = equations,
    latent = "L",
    latent_method = "explicit",
    n.iter = 2000,
    n.burnin = 1000,
    quiet = FALSE
)

cat("\nParameters monitored (Explicit):\n")
print(names(fit_explicit$samples[[1]][1, ]))

cat("\nSummary (should show beta_X_L, beta_Y_L):\n")
print(summary(fit_explicit))

cat("\n✓ Both methods completed successfully!\n")
