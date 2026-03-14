library(because)
library(ape)

set.seed(123)
N <- 30
tree <- rtree(N)

cat("=== Test 1: Auto-detection with standard _se naming ===\n")
data1 <- data.frame(
    BodyMass = rnorm(N, mean = 50, sd = 10),
    BodyMass_se = runif(N, 0.5, 2), # Standard naming
    BrainSize = rnorm(N, mean = 100, sd = 20)
)

fit1 <- because(
    data = data1,
    tree = tree,
    equations = list(BrainSize ~ BodyMass),
    n.chains = 2,
    n.iter = 500,
    n.burnin = 250,
    quiet = TRUE
)
cat("✓ Auto-detection worked!\n\n")


cat("=== Test 2: Custom column name with se_col ===\n")
data2 <- data.frame(
    BodyMass = rnorm(N, mean = 50, sd = 10),
    BodyMass_SD = runif(N, 0.5, 2), # Non-standard name!
    BrainSize = rnorm(N, mean = 100, sd = 20)
)

fit2 <- because(
    data = data2,
    tree = tree,
    equations = list(BrainSize ~ BodyMass),
    variability = list(
        BodyMass = list(type = "se", se_col = "BodyMass_SD") # Custom column!
    ),
    n.chains = 2,
    n.iter = 500,
    n.burnin = 250,
    quiet = TRUE
)
cat("✓ Custom se_col worked!\n\n")


cat("=== Test 3: Custom column name for repeated measures ===\n")
data3 <- list(
    BodyMass_measurements = matrix(
        rnorm(N * 5, mean = 50, sd = 10),
        nrow = N,
        ncol = 5
    ),
    BrainSize = rnorm(N, mean = 100, sd = 20)
)

fit3 <- because(
    data = data3,
    tree = tree,
    equations = list(BrainSize ~ BodyMass),
    variability = list(
        BodyMass = list(type = "reps", obs_col = "BodyMass_measurements") # Custom!
    ),
    n.chains = 2,
    n.iter = 500,
    n.burnin = 250,
    quiet = TRUE
)
cat("✓ Custom obs_col worked!\n\n")


cat("=== All tests passed! ===\n")
cat("Summary of fit2:\n")
print(summary(fit2))
