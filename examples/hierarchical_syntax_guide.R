library(because)

cat("=======================================================\n")
cat("   HIERARCHICAL DATA - LME4 STYLE SYNTAX\n")
cat("=======================================================\n\n")

# Create test data with COLUMN names that match grouped
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

cat("✅ RECOMMENDED: Explicit hierarchy\n")
cat("===================================\n\n")

result <- suppressMessages({
    because(
        data = list(
            individual = individual_data,
            site_year = site_year_data
        ),
        levels = levels_spec,
        hierarchy = "site_year > individual", # Explicit - clearest!
        link_vars = c("site", "year"),
        equations = equations,
        n.chains = 1,
        n.iter = 100,
        quiet = TRUE
    )
})

cat("SUCCESS!\n\n")

cat("=======================================================\n")
cat("IMPORTANT DISTINCTION:\n")
cat("======================================================\n\n")

cat("Two different concepts:\n\n")

cat("1. HIERARCHICAL DATA LEVELS (dataset names):\n")
cat("   hierarchy = \"site_year > individual\"\n")
cat("   Tells because which DATASET to use for each d-sep test\n\n")

cat("2. RANDOM EFFECTS (grouping variables / columns):\n")
cat("   random = ~(1|site) + (1|year)\n")
cat("   Tells because which COLUMNS define random grouping\n\n")

cat("These are separate specifications!\n\n")

cat("CORRECT usage:\n")
cat("  data = list(individual = ..., site_year = ...)  # Dataset names\n")
cat("  hierarchy = \"site_year > individual\"            # Dataset hierarchy\n")
cat("  random = ~(1|site) + (1|year)                   # Column names\n\n")

cat("NESTING in random effects refers to COLUMNS:\n")
cat("  random = ~(1|site/plot)  means plot nested in site\n")
cat("  This expands to: (1|site) + (1|site:plot)\n")
cat("  Both 'site' and 'plot' must be COLUMNS in your data!\n\n")

cat("=======================================================\n")
