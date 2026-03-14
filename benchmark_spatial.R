library(devtools)
load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
load_all(
    "/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.occupancy"
)
library(spOccupancy)

set.seed(42)

# --- 1. SIMULATE SPATIAL DATA ---
N_sites <- 100
J_visits <- 3

# Coordinates
coords <- cbind(runif(N_sites, 0, 10), runif(N_sites, 0, 10))
dist_mat <- as.matrix(dist(coords))

# Spatial Random Effect (Exponential Kernel)
phi <- 3 / 5 # effective range of 5 units
sigma_sq <- 1
Sigma <- sigma_sq * exp(-dist_mat * phi)
w <- MASS::mvrnorm(1, rep(0, N_sites), Sigma)

# Occupancy
beta0 <- 0.5
psi <- plogis(beta0 + w)
z <- rbinom(N_sites, 1, psi)

# Detection
p <- 0.6
y <- matrix(
    rbinom(N_sites * J_visits, 1, z * p),
    nrow = N_sites,
    ncol = J_visits
)

# Format for because
df <- data.frame(
    Site = sprintf("site_%03d", 1:N_sites),
    coord_x = coords[, 1],
    coord_y = coords[, 2]
)
df$Y <- I(y)

# Format for spOccupancy
data_sp <- list(
    y = y,
    occ.covs = data.frame(int = rep(1, N_sites)),
    det.covs = data.frame(int = rep(1, N_sites)),
    coords = coords
)

# --- 2. FIT JAGS (via because) ---
message("\nFitting JAGS (because)...")
time_jags <- system.time({
    fit_jags <- because(
        list(Y ~ 1, p_Y ~ 1),
        data = df,
        family = c(Y = "occupancy"),
        structure = Sigma,
        n.iter = 20000,
        n.burnin = 5000,
        n.chains = 3,
        parallel = TRUE,
        n.cores = 3,
        quiet = TRUE,
        DIC = FALSE,
        WAIC = FALSE
    )
})
message(sprintf("JAGS Elapsed: %.2f s", time_jags["elapsed"]))

# --- 3. FIT NIMBLE (via because, marginalized) ---
message("\nFitting NIMBLE (because, marginalized)...")
time_nimble <- system.time({
    fit_nimble <- because(
        list(Y ~ 1, p_Y ~ 1),
        data = df,
        family = c(Y = "occupancy"),
        structure = Sigma,
        engine = "nimble",
        n.iter = 20000,
        n.burnin = 5000,
        n.chains = 3,
        parallel = TRUE,
        n.cores = 3,
        quiet = TRUE,
        DIC = FALSE,
        WAIC = FALSE
    )
})
message(sprintf(
    "NIMBLE (Marginalized) Elapsed: %.2f s",
    time_nimble["elapsed"]
))

# --- 4. FIT spOccupancy (spOcc) ---
message("\nFitting spOccupancy (spPGOcc)...")
time_sp <- system.time({
    fit_sp <- spPGOcc(
        occ.formula = ~1,
        det.formula = ~1,
        data = data_sp,
        n.batch = 400,
        batch.length = 50, # total 20000 samples
        n.burn = 5000,
        n.chains = 3,
        n.report = 100,
        NNGP = TRUE,
        n.neighbors = 10,
        verbose = FALSE
    )
})
message(sprintf("spOccupancy Elapsed: %.2f s", time_sp["elapsed"]))

# --- 5. RESULTS COMPARISON ---
message("\n--- PARAMETER ESTIMATES (Intercept beta0) ---")
cat("True beta0:", beta0, "\n")

# Get results
jags_sams <- as.matrix(fit_jags$samples)[, "alpha_Y"]
nimble_sams <- as.matrix(fit_nimble$samples)[, "alpha_Y"]
sp_sams <- as.matrix(fit_sp$beta.samples)

cat("JAGS Mean:", mean(jags_sams), " (SD:", sd(jags_sams), ")\n")
cat("NIMBLE Mean:", mean(nimble_sams), " (SD:", sd(nimble_sams), ")\n")
cat("spOccupancy Mean:", mean(sp_sams), " (SD:", sd(sp_sams), ")\n")

message("\n--- CONVERGENCE (R-hat) ---")
cat(
    "JAGS R-hat (alpha_Y):",
    coda::gelman.diag(fit_jags$samples[, "alpha_Y"])$psrf[1],
    "\n"
)
cat(
    "NIMBLE R-hat (alpha_Y):",
    coda::gelman.diag(fit_nimble$samples[, "alpha_Y"])$psrf[1],
    "\n"
)

message("\n--- EFFICIENCY (samples/sec) ---")
message(sprintf("JAGS: %.2f", 60000 / time_jags["elapsed"]))
message(sprintf("NIMBLE: %.2f", 60000 / time_nimble["elapsed"]))
message(sprintf("spOccupancy: %.2f", 60000 / time_sp["elapsed"]))
