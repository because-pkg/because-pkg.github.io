library(because)

# Create example data matching user's format with different numbers of replicates
set.seed(123)
data_long <- data.frame(
    SP = c(rep("s1", 10), rep("s2", 7), rep("s3", 5)), # Different reps per species!
    BM = rnorm(22, mean = 0, sd = 1),
    NL = rnorm(22, mean = 0, sd = 1)
)

# Add some missing values
data_long$BM[c(5, 15)] <- NA
data_long$NL[c(8, 18)] <- NA

cat("\n=== Original Data ===\n")
cat("Species counts:\n")
print(table(data_long$SP))
print(head(data_long, 12))

# Create a simple tree
tree <- ape::read.tree(text = "(s1:1,s2:1,s3:1);")

cat("\n\n=== Testing because_format_data ===\n")
data_list <- because_format_data(data_long, species_col = "SP", tree = tree)

cat("\nFormatted data structure:\n")
str(data_list)

cat("\n\n=== BM Matrix ===\n")
print(data_list$BM)

cat("\n\n=== NL Matrix ===\n")
print(data_list$NL)

cat("\n\n=== Testing with because ===\n")
fit <- because(
    data = data_list,
    tree = tree,
    equations = list(NL ~ BM),
    n.iter = 1000,
    n.burnin = 500,
    quiet = TRUE
)

cat("\n✓ Model fitted successfully!\n")
print(summary(fit))
