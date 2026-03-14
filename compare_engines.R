library(because)
library(ape)
devtools::load_all()

set.seed(123)
n_species <- 50
tree <- rtree(n_species)
tree$tip.label <- paste0("sp", 1:n_species)

# Simulate phylogenetic signal
vcv_mat <- vcv(tree)
vcv_mat <- vcv_mat / max(vcv_mat) # standardize

# Cholesky decomposition for simulation
L <- t(chol(vcv_mat))
u <- L %*% rnorm(n_species)

# Predictor
X <- rnorm(n_species)
# Response with phylogenetic signal + predictor effect + noise
Y <- 2 + 1.5 * X + as.numeric(u) * 2 + rnorm(n_species, 0, 0.5)

dat <- data.frame(
    species = tree$tip.label,
    Y = Y,
    X = X
)

eqs <- list(Y ~ X)

message("\n======================================")
message("Fitting model with JAGS (Default)...")
message("======================================")
time_jags <- system.time({
    fit_jags <- because(
        equations = eqs,
        data = dat,
        tree = tree,
        engine = 'jags',
        n.iter = 500,
        n.burnin = 100,
        n.chains = 2,
        quiet = TRUE
    )
})

message("\n======================================")
message("Fitting model with NIMBLE...")
message("======================================")
time_nimble <- system.time({
    fit_nimble <- because(
        equations = eqs,
        data = dat,
        tree = tree,
        engine = 'nimble',
        n.iter = 500,
        n.burnin = 100,
        n.chains = 2,
        quiet = TRUE
    )
})

message("\n======================================")
message("RESULTS COMPARISON")
message("======================================")

cat("\n--- PARAMETER ESTIMATES (Mean) ---\n")
means_df <- data.frame(
    JAGS = fit_jags$summary$statistics[, "Mean"],
    NIMBLE = fit_nimble$summary$statistics[, "Mean"]
)
print(means_df)

cat("\n--- PARAMETER ESTIMATES (SD) ---\n")
sd_df <- data.frame(
    JAGS = fit_jags$summary$statistics[, "SD"],
    NIMBLE = fit_nimble$summary$statistics[, "SD"]
)
print(sd_df)

cat("\n--- EXECUTION TIME ---\n")
time_df <- data.frame(
    JAGS_seconds = time_jags["elapsed"],
    NIMBLE_seconds = time_nimble["elapsed"]
)
rownames(time_df) <- "Total Elapsed"
print(time_df)
