# Test Hierarchical Data Support (N Inflation Prevention)
# This tests whether separate data tables are correctly preserved
# and passed to JAGS with appropriate sample sizes.

devtools::load_all(".")
library(ape)

set.seed(123)

# --- 1. Simulate Hierarchical Data ---

# Level 1: Species (Top, coarsest)
N_Species <- 10
tree <- rtree(N_Species)
species_names <- tree$tip.label
tree$node.label <- NULL

# Species-level traits
Traits_df <- data.frame(
    Species = species_names,
    Trait = rnorm(N_Species) # e.g., Body Size
)

# Level 2: Observations (Finest)
# Each species observed multiple times
N_Obs <- 50
obs_species <- sample(species_names, N_Obs, replace = TRUE)
Y_df <- data.frame(
    Species = obs_species,
    Y = rnorm(N_Obs) # Simple Gaussian response for testing
)

# --- 2. Test Auto-Detection (SIMPLIFIED SYNTAX) ---

print("Testing auto-detect hierarchical structure...")

hierarchical_data <- list(
    Traits = Traits_df,
    Y = Y_df
)

# Get all equation variables
eq_vars <- c("Trait", "Y")

# Call auto-detection
auto_result <- auto_detect_hierarchical(
    hierarchical_data,
    eq_vars,
    quiet = FALSE
)

print("Auto-detected result:")
print(auto_result)

# Verify auto-detection
stopifnot(
    "Trait should be in Traits level" = "Trait" %in% auto_result$levels$Traits
)
stopifnot("Y should be in Y level" = "Y" %in% auto_result$levels$Y)
stopifnot(
    "Species should be detected as link var" = "Species" %in%
        auto_result$link_vars
)
stopifnot(
    "Hierarchy should show Traits > Y" = grepl(
        "Traits.*>.*Y",
        auto_result$hierarchy
    )
)

print("SUCCESS: Auto-detection works correctly!")

# --- 3. Test Data Preparation ---

print("Testing prepare_hierarchical_jags_data...")

prep_result <- prepare_hierarchical_jags_data(
    hierarchical_info = list(
        data = hierarchical_data,
        levels = auto_result$levels,
        hierarchy = auto_result$hierarchy,
        link_vars = auto_result$link_vars
    ),
    vars_needed = eq_vars
)

print("Sample sizes extracted:")
print(prep_result$n_vec)

# Verify N_Traits is NOT inflated
stopifnot("N_Traits should be 10" = prep_result$n_vec$N_Traits == N_Species)
stopifnot("N_Y should be 50" = prep_result$n_vec$N_Y == N_Obs)

print("SUCCESS: Sample sizes are correct (no N inflation).")

# --- 4. Test Model String Generation ---

print("Testing model string generation...")

model_output <- because_model(
    equations = list(Y ~ Trait),
    family = NULL,
    hierarchical_info = list(
        data = hierarchical_data,
        levels = auto_result$levels,
        hierarchy = auto_result$hierarchy,
        link_vars = auto_result$link_vars
    )
)

# Check for hierarchical loop
has_n_y_loop <- grepl("for \\(i in 1:N_Y\\)", model_output$model)
has_trait_indexing <- grepl(
    "Trait\\[Traits_idx_Y\\[i\\]\\]",
    model_output$model
)

stopifnot("Model should use N_Y loop" = has_n_y_loop)
stopifnot("Model should use cross-level indexing" = has_trait_indexing)

print("SUCCESS: Model string uses hierarchical bounds and indexing!")

print("=== All hierarchical data tests passed! ===")
