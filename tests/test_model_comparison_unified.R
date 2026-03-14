library(because)
library(ape)

# Create synthetic data
set.seed(42)
N <- 40
tree <- rtree(N)
df <- data.frame(
    SP = tree$tip.label,
    Y = rnorm(N, 10, 2),
    X1 = rnorm(N),
    X2 = rnorm(N)
)

cat("\n=== TEST 1: Run and Compare Specs (Mode 2) ===\n")
specs <- list(
    m1 = list(equations = list(Y ~ X1)),
    m2 = list(equations = list(Y ~ X1 + X2))
)

# Test positional arguments for specs, data, tree
res_batch <- because_compare(specs, df, tree, n.cores = 1, n.iter = 1000)
print(res_batch$comparison)

cat("\n=== TEST 2: Compare Fitted Models (Mode 1) ===\n")
fit1 <- res_batch$results$m1
fit2 <- res_batch$results$m2

# Test positional arguments for fits
comp_fits <- because_compare(fit1, fit2)
print(comp_fits)

cat("\n✓ SUCCESS\n")
