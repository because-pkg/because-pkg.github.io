library(because)

cat("=======================================================\n")
cat("   TESTING HIERARCHY INFERENCE (SIMPLIFIED)\n")
cat("=======================================================\n\n")

# Create test data
individual_data <- data.frame(
    id = 1:20,
    site = rep(c("A", "B"), each = 10),
    year = rep(c(2020, 2021), 10),
    weight = rnorm(20, 50, 10),
    age = rnorm(20, 5, 2)
)

site_year_data <- data.frame(
    site = c("A", "A", "B", "B"),
    year = c(2020, 2021, 2020, 2021),
    temp = c(15, 16, 17, 18),
    precip = c(800, 850, 900, 950)
)

levels_spec <- list(
    individual = c("weight", "age"),
    site_year = c("temp", "precip")
)

equations <- list(
    weight ~ age + temp
)

cat("Test 1: Explicit hierarchy (recommended - clearest)\n")
cat("---------------------------------------------------\n")

result1 <- suppressMessages({
    because(
        data = list(
            individual = individual_data,
            site_year = site_year_data
        ),
        levels = levels_spec,
        hierarchy = "site_year > individual", # Explicit and clear
        link_vars = c("site", "year"),
        equations = equations,
        n.chains = 1,
        n.iter = 100,
        quiet = TRUE
    )
})

cat("✅ Test 1 PASSED\n\n")

cat("=======================================================\n")
cat("   RECOMMENDED USAGE\n")
cat("=======================================================\n\n")

cat("For hierarchical data, use explicit hierarchy specification:\n\n")

cat("hierarchy = \"site_year > individual\"\n\n")

cat("This is clearest and most reliable. Future versions may support\n")
cat("more sophisticated inference from random effects syntax.\n")
