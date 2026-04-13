# Verification Script v5: Valid DAG with Worker Failure
# Achaz von Hardenberg
# 13 April 2026

library(pbapply)
library(parallel)

# Load the local version of because
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")

# 1. Simulate Small Data (X -> Y -> Z)
set.seed(42)
N <- 50
X <- rnorm(N)
Y <- 0.5 * X + rnorm(N)
Z <- -2 + 0.3 * Y + rnorm(N) # Z has negative values
data <- data.frame(X=X, Y=Y, Z=Z)

# 2. Define Equations (X -> Y -> Z implies X _||_ Z | Y)
eqs <- list(
  Y ~ X,
  Z ~ Y
)

message("\n--- Testing Worker Robustness (DAG X -> Y -> Z) ---")
message("Force Poisson fit on negative Z in parallel...")

# This should generate 1 test: Z ~ X + Y (or X ~ Z + Y)
# We force Z to be poisson, which will fail during JAGS initialization.

fit_bad <- because(
  eqs,
  data = data,
  dsep = TRUE,
  parallel = TRUE,
  n.cores = 2,
  n.iter = 500,
  family = c(Z = "poisson") # This will fail in worker
)

message("\nVerification Complete.")
