library(because)
library(coda)

# Simulate data similar to rhino structure
set.seed(42)
N <- 50

# Create correlated structure from latent variable
L1 <- rnorm(N)
LS <- 0.6 * L1 + rnorm(N, 0, 0.3)
NL <- 0.5 * L1 + 0.4 * rnorm(N) + rnorm(N, 0, 0.3) # RS effect + L1 effect
RS <- rnorm(N)
NL <- NL + 0.4 * RS # Add RS effect
DD <- 0.7 * NL + rnorm(N, 0, 0.3)

# Standardize for better convergence
df <- data.frame(
    SP = paste0("sp", 1:N),
    LS = scale(LS)[, 1],
    NL = scale(NL)[, 1],
    RS = scale(RS)[, 1],
    DD = scale(DD)[, 1]
)

equations <- list(
    LS ~ L1,
    NL ~ L1 + RS,
    DD ~ NL
)

cat("Running MAG model with induced correlations...\n")

fit <- because(
    data = df,
    id_col = "SP",
    latent = "L1",
    equations = equations,
    latent_method = "correlations",
    n.chains = 3,
    n.iter = 3000,
    n.burnin = 1000,
    quiet = TRUE
)

cat("\n=== JAGS Model Check ===\n")
model_str <- fit$model_code
cat("Looking for 'tau_obs_NL ~ dgamma' (should be estimated):\n")
if (grepl("tau_obs_NL ~ dgamma", model_str)) {
    cat("✓ PASS: tau_obs is estimated (not fixed)\n")
} else if (grepl("tau_obs_NL <- 10000", model_str)) {
    cat("✗ FAIL: tau_obs is still fixed at 10000\n")
}

cat("\nLooking for 'RS[i] ~ dnorm' (should NOT be present):\n")
if (grepl("RS\\[i\\] ~ dnorm", model_str)) {
    cat("✗ FAIL: RS is being estimated as stochastic\n")
} else {
    cat("✓ PASS: RS is treated as observed data\n")
}

cat("\n=== Convergence Diagnostics ===\n")
mcmc_samples <- fit$samples
rhats <- gelman.diag(mcmc_samples)$psrf[, 1]
beta_rhat <- rhats["beta_NL_RS"]
cat(sprintf(
    "beta_NL_RS Rhat: %.3f %s\n",
    beta_rhat,
    ifelse(beta_rhat < 1.1, "✓ GOOD", "✗ POOR")
))

neff <- effectiveSize(mcmc_samples)["beta_NL_RS"]
cat(sprintf(
    "beta_NL_RS n.eff: %.1f %s\n",
    neff,
    ifelse(neff > 100, "✓ GOOD", ifelse(neff > 20, "~ OK", "✗ POOR"))
))

cat("\n=== Summary ===\n")
print(summary(mcmc_samples[, "beta_NL_RS"]))
