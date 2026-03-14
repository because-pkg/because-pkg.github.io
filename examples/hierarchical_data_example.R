# Hierarchical Data Example - Current Functionality
# This demonstrates the hierarchical data structure feature (in progress)

library(because)

# ============================================================================
# CREATE EXAMPLE DATA
# ============================================================================

# Individual-level data (20 marmots across 2 sites, 2 years)
individual_data <- data.frame(
    individual_id = 1:20,
    site = rep(c("Valley", "Ridge"), each = 10),
    year = rep(c(2020, 2021), 10),
    weight = rnorm(20, 3500, 500), # grams
    age = sample(1:8, 20, replace = TRUE),
    activity = rnorm(20, 50, 10) # minutes active per day
)

# Site-year level data (environmental conditions - only 4 unique values!)
site_year_data <- data.frame(
    site = c("Valley", "Valley", "Ridge", "Ridge"),
    year = c(2020, 2021, 2020, 2021),
    mean_temp = c(15.2, 16.1, 12.8, 13.5), # °C
    total_precip = c(820, 780, 950, 1020), # mm
    growing_season = c(120, 125, 110, 115) # days
)

cat("Individual-level data:\n")
print(head(individual_data, 3))
cat("\nSite-year level data (only 4 rows for 20 individuals!):\n")
print(site_year_data)

# ============================================================================
# DEFINE HIERARCHICAL STRUCTURE
# ============================================================================

# Specify which variables belong to which level
levels_spec <- list(
    individual = c("weight", "age", "activity"),
    site_year = c("mean_temp", "total_precip", "growing_season")
)

# Define causal model
equations <- list(
    weight ~ age + mean_temp, # Individual outcome ~ individual + site-year predictors
    activity ~ age + growing_season # Another mixed-level relationship
)

# ============================================================================
# RUN MODEL WITH HIERARCHICAL DATA
# ============================================================================

cat("\n=== Running because with Hierarchical Data ===\n\n")

result <- because(
    data = list(
        individual = individual_data,
        site_year = site_year_data
    ),
    levels = levels_spec,
    hierarchy = "site_year > individual",
    link_vars = c("site", "year"),
    equations = equations,
    n.chains = 2,
    n.iter = 1000,
    quiet = FALSE
)

cat("\n=== Model Results ===\n")
summary(result)

# ============================================================================
# WHAT'S WORKING
# ============================================================================

cat("\n=== ✅ WHAT'S WORKING ===\n\n")
cat("1. Data Validation:\n")
cat("   - Checks that all datasets are data.frames ✓\n")
cat("   - Verifies variables exist in specified levels ✓\n")
cat("   - Validates hierarchy specification ✓\n")
cat("   - Ensures link variables exist in all datasets ✓\n\n")

cat("2. Data Assembly:\n")
cat("   - Automatically joins site_year and individual data ✓\n")
cat("   - Uses appropriate linking variables (site, year) ✓\n")
cat("   - Creates full dataset for model estimation ✓\n")
cat("   - Sample size:", nrow(individual_data), "observations ✓\n\n")

cat("3. Model Estimation:\n")
cat("   - Mixed-level predictors work correctly ✓\n")
cat("   - Environmental variables (mean_temp, etc.) properly joined ✓\n")
cat("   - MCMC sampling completes successfully ✓\n\n")

# ============================================================================
# CURRENT LIMITATIONS
# ============================================================================

cat("\n=== ⚠️  CURRENT LIMITATIONS ===\n\n")
cat("D-separation tests not yet optimized for hierarchical data:\n")
cat("  - Tests between ONLY site-year variables still use all 20 rows\n")
cat("  - Should use aggregated data (4 rows) for those tests\n")
cat("  - This causes sample size inflation for higher-level relationships\n\n")

cat("Example: Testing 'mean_temp _||_ total_precip | growing_season'\n")
cat("  Current: Uses 20 observations (incorrect - inflated)\n")
cat(
    "  Should use: 4 observations (correct - actual site-year combinations)\n\n"
)

# ============================================================================
# WORKAROUND FOR D-SEPARATION
# ============================================================================

cat("=== WORKAROUND (until Phase 3 complete) ===\n\n")
cat("For d-separation tests between site-year variables only:\n")
cat("Manually use aggregated data:\n\n")

cat("# Test relationship between environmental variables\n")
cat("fit_climate <- because(\n")
cat("  equations = list(total_precip ~ mean_temp),\n")
cat("  data = site_year_data,  # Use aggregated data directly\n")
cat("  dsep = TRUE\n")
cat(")\n\n")

# ============================================================================
# NEXT STEPS
# ============================================================================

cat("=== NEXT IMPLEMENTATION STEPS ===\n\n")
cat("Phase 3 (in progress):\n")
cat("  - Analyze variables in each d-sep test\n")
cat("  - Select appropriate dataset (aggregated vs full)\n")
cat("  - Pass correct data to each test\n")
cat("  - Maintain proper sample sizes\n\n")

cat("Once complete, d-separation will automatically:\n")
cat("  ✓ Use site_year_data for tests between environmental variables\n")
cat("  ✓ Use full data for tests involving individual-level variables\n")
cat("  ✓ Prevent sample size inflation\n")
cat("  ✓ Produce statistically correct conditional independence tests\n")
