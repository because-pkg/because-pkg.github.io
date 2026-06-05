lines <- readLines("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because/R/because.R")

# Find the lines where family is passed to because_py$fit
idx1 <- grep("family = family,", lines)

for (i in idx1) {
  lines[i] <- sub("family = family,", "family = if (!is.null(family)) as.list(family) else NULL,", lines[i])
}

writeLines(lines, "/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because/R/because.R")
