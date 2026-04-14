# Thorson (2025) Benchmark: OU Comparative Analysis
# This script executes the benchmark using the 'because' package and the custom OU extension.

# 1. Environment Setup
library(because)
source("scratch/thorson_benchmark/ou_extension.R")

# 2. Load Harmonized Data
cat("Loading cleaned data...\n")
data_list <- readRDS("data/thorson_benchmark/cleaned_data.rds")
traits <- data_list$traits
tree <- data_list$tree

# 3. Define the OU Selection Suite (Alpha values to test via BMA)
alphas <- c(0.001, 0.01, 0.1, 1.0, 10.0) 

# Prepare OU structures for each trait
cat("Preparing OU precision matrices for BMA...\n")
# mass_ou, metabolism_ou, range_ou
mass_structure       <- ou_structure(tree, trait_name = "log_mass", alphas = alphas)
metabolism_structure <- ou_structure(tree, trait_name = "log_metabolism", alphas = alphas)
range_structure      <- ou_structure(tree, trait_name = "log_range", alphas = alphas)

# 4. Define the Causal Model (GGMM)
# We model the three traits with their respective phylogenetic OU structures.
# Thorson (2025) highlights that range has higher alpha than mass.
# Note: In because(), we use the 'equations' argument (list or string).
model_equations <- "log_mass ~ 1
                    log_metabolism ~ log_mass
                    log_range ~ log_mass"

# Execute because() with BMA across OU kernels
cat("Running Bayesian OU Path Analysis (BMA)...\n")
# Note: For speed in this demonstration, we use a subset of the data.
# 500 tips with 3D precision matrices in JAGS is perfectly feasible.
set.seed(42)
subset_idx <- sample(1:nrow(traits), 500)
traits_sub <- traits[subset_idx, ]
tree_sub <- ape::keep.tip(tree, traits_sub$species)

# Re-prepare structures for the subset
# We name the list elements to match the response variables
# This prefixes the K index in JAGS (e.g. K_log_mass)
cat("Computing OU precision matrices for subset...\n")
mass_struct_sub  <- ou_structure(tree_sub, "log_mass", alphas)
range_struct_sub <- ou_structure(tree_sub, "log_range", alphas)

result <- because(
  equations = model_equations,
  data = traits_sub,
  structure = list(
    log_mass  = mass_struct_sub,
    log_range = range_struct_sub
  ),
  n.iter = 5000,
  n.burnin = 2000,
  quiet = FALSE,
  verbose = TRUE
)

# 5. Extract Selection Strength (BMA Posterior Probabilities)
# We look for the K_... variables in the samples
cat("\nExtracting BMA results...\n")

get_bma_probs <- function(result, var_name, alphas) {
  # Extract samples for the selection index K
  k_samples <- result$samples[, paste0("K_", var_name), drop=FALSE]
  k_samples <- as.vector(as.matrix(k_samples))
  
  # Frequency count
  probs <- table(factor(k_samples, levels = seq_along(alphas))) / length(k_samples)
  names(probs) <- as.character(alphas)
  return(probs)
}

post_alpha_mass  <- get_bma_probs(result, "log_mass", alphas)
post_alpha_range <- get_bma_probs(result, "log_range", alphas)

cat("\n--- Benchmark Results ---\n")
cat("Alpha Suite:", paste(alphas, collapse=", "), "\n")
cat("\nPosterior Probs for Body Mass Selection (Alpha):\n")
print(round(post_alpha_mass, 3))
cat("\nPosterior Probs for Geographic Range Selection (Alpha):\n")
print(round(post_alpha_range, 3))

# Check which alpha is the 'Mode' (Highest Posterior Probability)
map_mass <- alphas[which.max(post_alpha_mass)]
map_range <- alphas[which.max(post_alpha_range)]

cat(sprintf("\nConclusion: Body Mass follows alpha=%.3f, while Geographic Range follows alpha=%.3f\n", 
           map_mass, map_range))

# Save result for manuscript
saveRDS(result, "data/thorson_benchmark/benchmark_result.rds")
cat("\nAnalysis complete. Results saved to data/thorson_benchmark/benchmark_result.rds\n")
