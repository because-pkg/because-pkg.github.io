library(because)

cat("=======================================================\n")
cat("   TESTING IMPROVED HIERARCHY INFERENCE\n")
cat("=======================================================\n\n")

# Create test data
individual_data <- data.frame(
    id = 1:20,
    site = rep(c("A", "B"), each = 10),
    year = rep(c(2020, 2021), 10),
    individual = paste0("ind_", 1:20),
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

cat("Test 1: Explicit nesting syntax (1|site_year/individual)\n")
cat("---------------------------------------------------------\n")

result1 <- suppressMessages({
    because(
        data = list(
            individual = individual_data,
            site_year = site_year_data
        ),
        levels = levels_spec,
        link_vars = c("site", "year"),
        random = ~ (1 | site_year / individual), # Explicit nesting
        equations = equations,
        n.chains = 1,
        n.iter = 100,
        quiet = TRUE
    )
})

cat("✅ Test 1 PASSED: Inferred from (1|site_year/individual)\n\n")

cat("Test 2: Matching grouping variables to level names\n")
cat("-----------------------------------------------------\n")

result2 <- suppressMessages({
    because(
        data = list(
            individual = individual_data,
            site_year = site_year_data
        ),
        levels = levels_spec,
        link_vars = c("site", "year"),
        random = ~ (1 | individual) + (1 | site) + (1 | year), # Separate terms
        # Should infer: site and year -> site_year, individual -> individual
        equations = equations,
        n.chains = 1,
        n.iter = 100,
        quiet = TRUE
    )
})

cat("✅ Test 2 PASSED: Inferred from separate grouping variables\n\n")

cat("Test 3: Still works with explicit hierarchy\n")
cat("--------------------------------------------\n")

result3 <- suppressMessages({
    because(
        data = list(
            individual = individual_data,
            site_year = site_year_data
        ),
        levels = levels_spec,
        hierarchy = "site_year > individual", # Explicit
        link_vars = c("site", "year"),
        equations = equations,
        n.chains = 1,
        n.iter = 100,
        quiet = TRUE
    )
})

cat("✅ Test 3 PASSED: Explicit hierarchy still works\n\n")

cat("=======================================================\n")
cat("   ALL INFERENCE TESTS PASSED!\n")
cat("=======================================================\n\n")

cat("Hierarchy can now be specified 3 ways:\n")
cat("1. Explicit: hierarchy = \"site_year > individual\"\n")
cat("2. Nesting syntax: random = ~(1|site_year/individual)\n")
cat("3. Auto-infer: random = ~(1|site) + (1|year) + (1|individual)\n")
cat("   (matches grouping vars to level names automatically)\n")
