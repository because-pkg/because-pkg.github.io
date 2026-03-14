library(because)
library(ape)

# Test auto-detection of variability
set.seed(123)
N <- 30
tree <- rtree(N)

# Create data with SE (standard errors)
data <- data.frame(
    X = rnorm(N, mean = 5, sd = 2),
    X_se = runif(N, 0.1, 0.5), # Auto-detect this!
    Y = rnorm(N)
)
data$Y <- 0.5 * data$X + rnorm(N, sd = 0.5)

cat("=== Testing Auto-Detection of Standard Errors ===\n")
cat("Data columns:", paste(names(data), collapse = ", "), "\n\n")

# Run WITHOUT specifying variability argument
# Should auto-detect X_se
fit <- because(
    data = data,
    tree = tree,
    equations = list(Y ~ X),
    n.chains = 2,
    n.iter = 1000,
    n.burnin = 500,
    quiet = FALSE # Show messages to see auto-detection
)

cat("\n\n=== Model Ran Successfully ===\n")
cat("Parameters monitored:\n")
print(fit$monitor)

cat("\n\nSummary:\n")
print(summary(fit))
