# Test script for Generalized Covariance Support (Non-Phylogenetic & Custom Matrix)

# library(because)
devtools::load_all(".")
library(ape)
library(testthat)

# 1. Synthesize Independent Data
set.seed(123)
N <- 100
x1 <- rnorm(N)
beta_true <- 0.5
sigma_true <- 1.0
y <- beta_true * x1 + rnorm(N, 0, sigma_true)
data_df <- data.frame(y = y, x1 = x1)

# Equation
eq <- list(y ~ x1)

cat("--- Test 1: Independent Model (tree = NULL) ---\n")

# Run because with tree = NULL
fit_null <- because(
    data = data_df,
    tree = NULL,
    equations = eq,
    n.iter = 5000,
    n.burnin = 1000,
    n.thin = 5,
    quiet = FALSE
)

message("Model run complete.")

# Compare with simple linear regression (lm)
fit_lm <- lm(y ~ x1, data = data_df)
summ_lm <- summary(fit_lm)

# Extract estimates
beta_because <- mean(fit_null$samples[[1]][, "beta_y_x1"])
beta_lm <- coef(fit_lm)["x1"]

cat(sprintf("Beta (because): %.3f\n", beta_because))
cat(sprintf("Beta (lm):       %.3f\n", beta_lm))

diff <- abs(beta_because - beta_lm)
if (diff < 0.1) {
    cat("[PASS] Estimates match within tolerance.\n")
} else {
    cat("[FAIL] Estimates diverge significantly.\n")
}

cat("\n--- Test 2: Custom Matrix Support (Identity Matrix) ---\n")

# Pass diagonal matrix (Identity) - Should be equivalent to Null/Independent
# But with u_std generated (model interprets it as structured with Prec = Identity)
# Wait, if we pass Identity as Prec, u ~ N(0, I).
# Then Y = mu + u + epsilon.
# This is a model with two independent error terms (overparametrized if both are estimated freely).
# But let's just check if it RUNS.

identity_mat <- diag(N)
# Treat as Covariance Matrix
fit_matrix <- because(
    data = data_df,
    tree = identity_mat, # Pass matrix directly
    equations = eq,
    n.iter = 1000,
    quiet = TRUE
)

message("Matrix model run complete.")
if (!is.null(fit_matrix$samples)) {
    cat("[PASS] Custom matrix model ran successfully.\n")
} else {
    cat("[FAIL] Custom matrix model failed.\n")
}

cat("\n--- Test 3: Multinomial Independent (tree = NULL) ---\n")
# Brief test for multinomial syntax
N_cat <- 50
cats <- sample(1:3, N_cat, replace = TRUE)
x_cat <- rnorm(N_cat)
data_cat <- data.frame(cat = as.factor(cats), x = x_cat)

fit_multinom <- because(
    data = data_cat,
    tree = NULL,
    equations = list(cat ~ x),
    family = c(cat = "multinomial"),
    n.iter = 1000,
    quiet = FALSE
)

if (!is.null(fit_multinom$samples)) {
    cat("[PASS] Independent Multinomial model ran successfully.\n")
} else {
    cat("[FAIL] Independent Multinomial model failed.\n")
}
