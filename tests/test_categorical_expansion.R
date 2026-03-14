library(because)

library(ape)

set.seed(123)

# 1. Create test data with categorical predictor
N <- 20
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

cat("Categorical vars attribute:\n")
print(attr(data_list, "categorical_vars"))

# 3. Define equation with categorical variable
equations <- list(Y ~ Diet + Z)

cat("\n\nOriginal equation:\n")
print(equations[[1]])

# 4. Run because WITHOUT d-sep to see if expansion works
cat("\nRunning because with categorical variable in equation...\n")
fit <- because(
    data = data_list,
    tree = tree,
    equations = equations,
    n.iter = 500,
    n.burnin = 100,
    n.chains = 1,
    quiet = FALSE # Show expansion messages
)

cat("\n\n=== Model Code ===\n")
cat(fit$model_code)
cat("\n==================\n")

# Check monitored parameters
cat("\n\nMonitored parameters:\n")
print(colnames(fit$samples[[1]]))
