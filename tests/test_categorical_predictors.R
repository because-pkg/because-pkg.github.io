library(because)

library(ape)

set.seed(999)

# 1. Create test data with categorical predictor
N <- 50
tree <- rtree(N)

# Create categorical variable (3 levels)
diet_types <- c("Carnivore", "Herbivore", "Omnivore")
data_long <- data.frame(
    SP = tree$tip.label,
    BodyMass = rnorm(N, mean = 10, sd = 2),
    Diet = sample(diet_types, N, replace = TRUE), # Character variable
    Habitat = factor(sample(c("Forest", "Grassland"), N, replace = TRUE)) # Factor
)

cat("Original data (first 6 rows):\n")
print(head(data_long))
cat("\n")

# 2. Format data (should auto-create dummies)
cat("Formatting data with categorical variables...\n")
data_list <- because_format_data(data_long, species_col = "SP", tree = tree)

cat("\nFormatted data variables:\n")
print(names(data_list))
cat("\n")

# Check categorical metadata
if (!is.null(attr(data_list, "categorical_vars"))) {
    cat("Categorical variable info:\n")
    print(attr(data_list, "categorical_vars"))
    cat("\n")
}

# 3. Test with because using auto-expanding equations
cat("Running model with AUTO-EXPANDING categorical predictors...\n")
cat("Equation: BodyMass ~ Diet + Habitat (will auto-expand to dummies)\n\n")

fit <- because(
    data = data_list,
    tree = tree,
    equations = list(BodyMass ~ Diet + Habitat), # Auto-expands!
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 2,
    quiet = FALSE # Show expansion messages
)

cat("Model summary:\n")
summary_stats <- fit$summary$statistics
print(summary_stats[grep("beta", rownames(summary_stats)), ])

cat("\n✓ Categorical predictor test completed successfully!\n")
cat("\nInterpretation:\n")
cat(
    "- beta_BodyMass_Diet_Herbivore: Effect of Herbivore vs. Carnivore (reference)\n"
)
cat(
    "- beta_BodyMass_Diet_Omnivore: Effect of Omnivore vs. Carnivore (reference)\n"
)
cat(
    "- beta_BodyMass_Habitat_Grassland: Effect of Grassland vs. Forest (reference)\n"
)
