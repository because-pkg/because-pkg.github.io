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

cat("Running Model 1 (Simple)...\n")
fit1 <- because(
    data = df,
    structure = tree,
    id_col = "SP",
    equations = list(Y ~ X1),
    WAIC = TRUE,
    quiet = TRUE
)

cat("Running Model 2 (Complex)...\n")
fit2 <- because(
    data = df,
    structure = tree,
    id_col = "SP",
    equations = list(Y ~ X1 + X2),
    WAIC = TRUE,
    quiet = TRUE
)

cat("Running Model 3 (Null)...\n")
fit3 <- because(
    data = df,
    structure = tree,
    id_col = "SP",
    equations = list(Y ~ 1),
    WAIC = TRUE,
    quiet = TRUE
)

cat("\n=== Testing because_compare (individual args) ===\n")
comp <- because_compare(fit1, fit2, fit3)
print(comp)

cat("\n=== Testing because_compare (named args) ===\n")
comp_named <- because_compare(Simple = fit1, Complex = fit2, Null = fit3)
print(comp_named)

cat("\n=== Testing because_compare (list arg) ===\n")
comp_list <- because_compare(models = list(M1 = fit1, M2 = fit2, M3 = fit3))
print(comp_list)

cat("\n✓ SUCCESS!\n")
