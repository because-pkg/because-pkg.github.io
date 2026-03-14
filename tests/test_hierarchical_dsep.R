library(because)

# Test hierarchical data with d-separation
# This should show different sample sizes for different tests

# Individual level
individual_data <- data.frame(
    id = 1:20,
    site = rep(c("A", "B"), each = 10),
    year = rep(c(2020, 2021), 10),
    weight = rnorm(20, 50, 10),
    age = rnorm(20, 5, 2)
)

# Site-year level (only 4 unique combinations!)
site_year_data <- data.frame(
    site = c("A", "A", "B", "B"),
    year = c(2020, 2021, 2020, 2021),
    temp = c(15, 16, 17, 18),
    precip = c(800, 850, 900, 950)
)

# Equations with mixed levels
equations <- list(
    weight ~ age + temp, # Individual response ~ individual + site-year
    age ~ precip # Individual response ~ site-year
    # Implies potential test: weight _||_ precip | age, temp
    # Also: temp _||_ precip (both site-year - should use 4 obs!)
)

levels_spec <- list(
    individual = c("weight", "age"),
    site_year = c("temp", "precip")
)

cat("\n=== Testing D-separation with Hierarchical Data ===\n\n")

result <- because(
    data = list(
        individual = individual_data,
        site_year = site_year_data
    ),
    levels = levels_spec,
    hierarchy = "site_year > individual",
    link_vars = c("site", "year"),
    equations = equations,
    dsep = TRUE, # Run d-separation tests!
    n.chains = 2,
    n.iter = 500,
    quiet = FALSE
)

cat("\n=== SUCCESS! ===\n")
cat("If you saw different observation counts for different tests,\n")
cat("the hierarchical data selection is working!\n\n")

cat("Expected behavior:\n")
cat("- Tests with ONLY site-year variables (temp, precip): 4 observations\n")
cat("- Tests with individual variables: 20 observations\n")
cat("- Mixed tests: 20 observations (finest grain needed)\n")
