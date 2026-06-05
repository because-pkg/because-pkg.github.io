lines <- readLines("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because/R/because.R")

idx <- grep("Now family is a named character vector, family_objects stores originals", lines)
if (length(idx) > 0) {
  lines <- c(
    lines[1:(idx[1])],
    "    # Map 'bernoulli' to 'binomial' for internal consistency",
    "    family[family == \"bernoulli\"] <- \"binomial\"",
    lines[(idx[1]+1):length(lines)]
  )
  writeLines(lines, "/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because/R/because.R")
} else {
  print("Could not find insertion point!")
}
