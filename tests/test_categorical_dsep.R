library(because)

library(ape)

set.seed(123)

# 1. Create test data with categorical predictor
N <- 50
tree <- rtree(N)

# Create categorical variable (3 levels)
diet_types <- c("Carnivore", "Herbivore", "Omnivore")
data_long <- data.frame(
    SP = tree$tip.label,
    Y = rnorm(N),
    Z = rnorm(N),
    Diet = sample(diet_types, N, replace = TRUE)
)

# 2. Format data (auto-creates dummies: Diet_Herbivore, Diet_Omnivore)
data_list <- because_format_data(data_long, species_col = "SP", tree = tree)

# 3. Define model implying Y _||_ Diet | Z
# Model: Y ~ Z, Diet ~ Z
# Implied independence: Y _||_ Diet | Z
equations <- list(
    Y ~ Z,
    Diet ~ Z
)

cat("Running d-separation test with categorical variable 'Diet'...\n")
cat("Expected test: Y ~ Z + Diet\n")
cat(
    "Expected monitored params: betaY_Z, betaY_Diet_Herbivore, betaY_Diet_Omnivore\n\n"
)

# 4. Run because with dsep = TRUE
fit <- because(
    data = data_list,
    tree = tree,
    equations = equations,
    dsep = TRUE,
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 2,
    quiet = FALSE
)

# Print model code for debugging
cat("\n=== Generated JAGS Model Code ===\n")
print(fit$model_code)
cat("=================================\n")

# 5. Check results
dsep_results <- fit$dsep_tests
cat("\nNumber of d-sep tests:", length(dsep_results), "\n")

if (length(dsep_results) > 0) {
    cat("Test equation:", deparse(dsep_results[[1]]), "\n")

    # Check monitored parameters in the samples
    monitored_params <- colnames(fit$samples[[1]])
    cat("\nMonitored parameters:\n")
    print(monitored_params)

    # Verify dummy betas are present (allowing for suffix like _1)
    dummies_present <- c(
        any(grepl("betaDiet_Herbivore", monitored_params)),
        any(grepl("betaDiet_Omnivore", monitored_params))
    )

    if (all(dummies_present)) {
        cat("\n✓ SUCCESS: All categorical dummy coefficients monitored!\n")
    } else {
        cat("\n❌ FAILURE: Missing dummy coefficients in monitor list.\n")
    }
} else {
    cat("\n❌ FAILURE: No d-sep tests generated.\n")
}
