# Script to update because_model.R templates with robust replacements
setwd("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
lines <- readLines("R/because_model.R")

# 1. Define the property injection
injection <- c(
    "            # Structure properties",
    "            s_lvl <- get_struct_lvl(s_name, hierarchical_info)",
    "            s_bound <- if (is.null(s_lvl)) \"N\" else paste0(\"N_\", s_lvl)",
    "            s_zeros <- if (is.null(s_lvl)) \"zeros\" else paste0(\"zeros_\", s_lvl)"
)

# Find all occurrences of the structure loop start
new_lines <- c()
i <- 1
while (i <= length(lines)) {
    line <- lines[i]
    new_lines <- c(new_lines, line)

    if (grepl("for \\(s_name in structure_names\\) \\{", line)) {
        # Check if next line is already structure properties
        if (!grepl("Structure properties", lines[min(i + 1, length(lines))])) {
            new_lines <- c(new_lines, injection)
        }
    }
    i <- i + 1
}

# 2. Sequential replacements for the indices
content <- paste(new_lines, collapse = "\n")

# Replace [1:N, 1:N, K] -> [1:s_bound, 1:s_bound, K]
content <- gsub(
    "\\[1:N, 1:N, K\\]",
    "[1:\", s_bound, \", 1:\", s_bound, \", K]",
    content
)
content <- gsub(
    "\\[1:N, 1:N\\]",
    "[1:\", s_bound, \", 1:\", s_bound, \"]",
    content
)

# Replace zeros_obs etc with the structure zeros
# Note: some use get_zeros_name, some use get_var_level || "obs"
content <- gsub(
    "get_zeros_name\\(response, hierarchical_info\\)",
    "s_zeros",
    content
)
content <- gsub(
    "zeros_\"\\),\n\\s+get_var_level\\(response, hierarchical_info\\) %\\|\\| \"obs\"",
    "s_zeros, \"",
    content
)

# Replace the loop bound in dmnorm
content <- gsub(
    "\\[1:get_loop_bound\\(response, hierarchical_info\\)\\] ~ dmnorm",
    "[1:\", s_bound, \"] ~ dmnorm",
    content
)

# 3. Fix the mapping loop
content <- gsub(
    "\\[i\\] / sqrt\\(",
    "[\", get_struct_index(s_name, response, hierarchical_info), \"] / sqrt(",
    content
)

writeLines(content, "R/because_model.R")
cat("Successfully updated because_model.R\n")
