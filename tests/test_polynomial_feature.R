library(because)

cat("=== Testing Polynomial Transformation Feature ===\n\n")

# Simple test data
set.seed(123)
data <- data.frame(
    x = rnorm(50, 5, 2),
    y = rnorm(50, 10, 3)
)

# True relationship: y = 2 + 3*x + 0.5*x^2
data$y <- 2 + 3 * data$x + 0.5 * data$x^2 + rnorm(50, 0, 1)

cat("Test 1: Simple polynomial with I(x^2)\n")
cat("---------------------------------------\n")

fit1 <- because(
    equations = list(
        y ~ x + I(x^2) # Inline polynomial!
    ),
    data = data,
    n.chains = 1,
    n.iter = 500,
    quiet = FALSE
)

cat("\n\nJAGS Model Code (first 40 lines):\n")
cat(paste(head(strsplit(fit1$model, "\n")[[1]], 40), collapse = "\n"))

cat("\n\n=== Test Complete! ===\n")
cat("Check above for:\n")
cat("1. 'Detected 1 polynomial term(s): I(x^2)' message\n")
cat("2. JAGS code should have: x_pow2[i] <- x[i]^2\n")
cat("3. Model should estimate beta_x and beta_x_pow2\n")
