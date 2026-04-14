# Thorson (2025) Benchmark: Data Preparation
# This script cleans the PanTHERIA dataset and harmonizes it with the Mammal phylogeny.

library(ape)
library(dplyr)

# Paths (Local for portability)
trait_path <- "data/thorson_benchmark/PanTHERIA_1-0_WR05_Aug2008.txt"
tree_path  <- "data/thorson_benchmark/VertTree_mammals.tre"

# 1. Load Tree
cat("Loading Mammal Tree...\n")
tree <- read.tree(tree_path)

# 2. Load Trait Data
cat("Loading PanTHERIA Traits...\n")
# Note: PanTHERIA is tab-delimited and uses -999 for NA
traits_raw <- read.delim(trait_path, na.strings = "-999.00", stringsAsFactors = FALSE)

# 3. Clean and Select Core Traits (Thorson 2025 variables)
# - AdultBodyMass_g
# - BasalMetRate_mLO2hr
# - GR_Area_km2 (Geographic Range)
traits <- traits_raw %>%
  rename(
    species = MSW05_Binomial,
    mass = X5.1_AdultBodyMass_g,
    metabolism = X18.1_BasalMetRate_mLO2hr,
    range = X26.1_GR_Area_km2
  ) %>%
  select(species, mass, metabolism, range) %>%
  filter(!is.na(mass) | !is.na(metabolism) | !is.na(range))

# Standardize species names to match tree (Genus_species)
traits$species <- gsub(" ", "_", traits$species)

# 4. Harmonize Dataset and Tree
cat("Harmonizing species list...\n")
common_species <- intersect(tree$tip.label, traits$species)
cat("Found", length(common_species), "species present in both tree and traits.\n")

# Prune tree
pruned_tree <- keep.tip(tree, common_species)

# Filter and align traits
final_traits <- traits %>%
  filter(species %in% common_species) %>%
  arrange(match(species, pruned_tree$tip.label))

# 5. Log-transform (as per standard comparative methods)
final_traits <- final_traits %>%
  mutate(
    log_mass = log(mass),
    log_metabolism = log(metabolism),
    log_range = log(range)
  )

# 6. Save Clean Data
saveRDS(list(traits = final_traits, tree = pruned_tree), "data/thorson_benchmark/cleaned_data.rds")
cat("Cleaned data saved to data/thorson_benchmark/cleaned_data.rds\n")
