library(because)

# Simulate data instead of loading rhino.dat
set.seed(123)
N <- 50
rhino.tree <- ape::rtree(N)
rhino.tree$edge.length <- rhino.tree$edge.length /
    max(ape::branching.times(rhino.tree))

BM <- rnorm(N)
NL <- rnorm(N)
DD <- rnorm(N)
LS <- rnorm(N)
RS <- rnorm(N)

data_list <- list(
    BM = BM,
    NL = NL,
    DD = DD,
    LS = LS,
    RS = RS
)

equations_1 <- list(LS ~ BM, NL ~ BM, DD ~ NL, RS ~ DD)

# Test model saving
cat("\n=== Testing Model Saving for d-sep ===\n")
fit_rhino_dsep <- because(
    data = data_list,
    tree = rhino.tree,
    equations = equations_1,
    dsep = TRUE,
    n.iter = 500,
    n.burnin = 250,
    quiet = TRUE
)

cat("\n\n=== Checking Models Field ===\n")
cat("Number of models saved:", length(fit_rhino_dsep$models), "\n")
cat("Names of list elements:", names(fit_rhino_dsep), "\n\n")

cat("=== Model 1 (first 30 lines) ===\n")
model_1_lines <- strsplit(fit_rhino_dsep$models[[1]], "\n")[[1]]
cat(paste(model_1_lines[1:min(30, length(model_1_lines))], collapse = "\n"))

cat("\n\n=== Model 3 (first 30 lines) ===\n")
model_3_lines <- strsplit(fit_rhino_dsep$models[[3]], "\n")[[1]]
cat(paste(model_3_lines[1:min(30, length(model_3_lines))], collapse = "\n"))

cat("\n\nTest completed successfully!\n")
