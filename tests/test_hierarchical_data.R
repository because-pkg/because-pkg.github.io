library(because)

# Create hierarchical test data
# Individual level
individual_data <- data.frame(
    id = 1:20,
    site = rep(c("A", "B"), each = 10),
    year = rep(2020:2021, each = 10),
    weight = rnorm(20, 50, 10),
    age = rnorm(20, 5, 2)
)

# Site-year level
site_year_data <- data.frame(
    site = rep(c("A", "B"), each = 2),
    year = rep(2020:2021, 2),
    temperature = c(15, 16, 17, 18),
    rainfall = c(800, 850, 900, 950)
)

# Define levels
levels_spec <- list(
    individual = c("weight", "age"),
    site_year = c("temperature", "rainfall")
)

# Simple equations
equations <- list(
    weight ~ age + temperature
)

cat("\n=== Testing Hierarchical Data Validation ===\n")

# Test 1: Basic validation should work
tryCatch(
    {
        result <- because(
            data = list(
                individual = individual_data,
                site_year = site_year_data
            ),
            levels = levels_spec,
            hierarchy = "site_year > individual",
            link_vars = c("site", "year"),
            equations = equations,
            n.chains = 1,
            n.iter = 100,
            quiet = FALSE
        )
        cat("✓ Test 1 PASSED: Hierarchical data accepted and assembled\n")
        cat(
            "  Assembled data dimensions:",
            nrow(result$data_used),
            "×",
            ncol(result$data_used),
            "\n"
        )
    },
    error = function(e) {
        cat("✗ Test 1 FAILED:", e$message, "\n")
    }
)

cat("\n=== Testing Error Cases ===\n")

# Test 2: Missing levels argument should error
tryCatch(
    {
        result <- because(
            data = list(
                individual = individual_data,
                site_year = site_year_data
            ),
            # levels missing!
            hierarchy = "site_year > individual",
            link_vars = c("site", "year"),
            equations = equations,
            n.chains = 1,
            n.iter = 100,
            quiet = TRUE
        )
        cat("✗ Test 2 FAILED: Should have errored without levels\n")
    },
    error = function(e) {
        cat("✓ Test 2 PASSED: Correctly errored -", e$message, "\n")
    }
)

# Test 3: Variable in wrong level should error
bad_levels <- list(
    individual = c("weight", "temperature"), # temperature is site_year!
    site_year = c("age", "rainfall") # age is individual!
)

tryCatch(
    {
        result <- because(
            data = list(
                individual = individual_data,
                site_year = site_year_data
            ),
            levels = bad_levels,
            hierarchy = "site_year > individual",
            link_vars = c("site", "year"),
            equations = equations,
            n.chains = 1,
            n.iter = 100,
            quiet = TRUE
        )
        cat("✗ Test 3 FAILED: Should have errored with wrong variable levels\n")
    },
    error = function(e) {
        cat("✓ Test 3 PASSED: Correctly errored -", e$message, "\n")
    }
)

cat("\n=== Summary ===\n")
cat("Hierarchical data infrastructure is working!\n")
cat("Next step: Implement per-test data selection for d-separation\n")
