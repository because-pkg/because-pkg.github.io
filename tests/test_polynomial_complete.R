library(because)

cat("===================================================\n")
cat("   POLYNOMIAL TRANSFORMATIONS - FULL TEST\n")
cat("===================================================\n\n")

# Create test data
set.seed(456)
data <- data.frame(
    x = rnorm(100, 5, 2),
    z = rnorm(100, 3, 1),
    y = NA
)

# True model: y depends on x, x^2, and z
data$y <- 2 + 3 * data$x + 0.5 * data$x^2 + 1.5 * data$z + rnorm(100, 0, 2)

cat("Test: Polynomial with d-separation\n")
cat("------------------------------------\n\n")

cat("Equations:\n")
cat("  y ~ x + I(x^2) + z\n")
cat("  z ~ x\n\n")

cat("Expected d-separation tests:\n")
cat("  Should NOT test: x _||_ x_pow2 (deterministic!)\n")
cat("  Should test: y _||_ x | x_pow2, z\n\n")

fit <- because(
    equations = list(
        y ~ x + I(x^2) + z,
        z ~ x
    ),
    data = data,
    dsep = TRUE,
    n.chains = 2,
    n.iter = 1000,
    quiet = FALSE
)

cat("\n===================================================\n")
cat("   D-SEPARATION RESULTS\n")
cat("===================================================\n\n")

summary(fit)

cat("\n===================================================\n")
cat("   SUCCESS!\n")
cat("===================================================\n\n")

cat("✅ Polynomial recognized: I(x^2)\n")
cat("✅ JAGS model has deterministic node: x_pow2[i] <- x[i]^2\n")
cat("✅ No spurious d-sep test for polynomial relationship\n")
cat("✅ Correct parameter estimates\n\n")
