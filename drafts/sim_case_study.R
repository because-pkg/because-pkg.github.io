devtools::load_all(".")
devtools::load_all("../because.phybase")
library(ape)
library(MASS)

set.seed(42)

# ==========================================
# 1. Evolutionary Level: Species Traits
# ==========================================
N_species <- 50
tree <- rcoal(N_species)
tree$tip.label <- paste0("Sp_", 1:N_species)
tree$edge.length <- tree$edge.length / max(branching.times(tree))
vcv_mat <- vcv(tree)

lambda_true <- 0.7 # Phylogenetic signal for Thermal_Tol
phylo_effects <- mvrnorm(
    1,
    mu = rep(0, N_species),
    Sigma = lambda_true * vcv_mat
)
names(phylo_effects) <- tree$tip.label

Body_Mass <- rnorm(N_species, mean = 50, sd = 10)
names(Body_Mass) <- tree$tip.label

# Thermal Tolerance (Species-Level)
beta_tol_mass <- 0.8
Thermal_Tol <- beta_tol_mass *
    scale(Body_Mass)[, 1] +
    phylo_effects +
    rnorm(N_species, 0, 0.5)
names(Thermal_Tol) <- tree$tip.label

# ==========================================
# 2. Ecological Level: Site Characteristics
# ==========================================
N_sites <- 50
site_coords <- matrix(runif(N_sites * 2, 0, 100), ncol = 2)
rownames(site_coords) <- paste0("Site_", 1:N_sites)
dist_mat <- as.matrix(dist(site_coords))
V_spatial <- exp(-dist_mat / 20) # Spatial correlation matrix

site_effects <- mvrnorm(
    1,
    mu = rep(0, N_sites),
    Sigma = 0.5^2 * V_spatial
)
names(site_effects) <- rownames(site_coords)

Elevation <- rnorm(N_sites, mean = 1000, sd = 300)
names(Elevation) <- rownames(site_coords)

# Temperature (Site-Level)
beta_temp_elev <- -1.5
Temperature <- beta_temp_elev *
    scale(Elevation)[, 1] +
    site_effects +
    rnorm(N_sites, 0, 0.5)
names(Temperature) <- rownames(site_coords)

# ==========================================
# 3. Observation Level: Species Abundance
# ==========================================
N_obs <- 500
obs_df <- data.frame(
    phylo = sample(tree$tip.label, N_obs, replace = TRUE),
    Site = sample(rownames(site_coords), N_obs, replace = TRUE)
)

# Merge Hierarchical Data into Observation DataFrame
obs_df$Body_Mass <- Body_Mass[obs_df$phylo]
obs_df$Thermal_Tol <- Thermal_Tol[obs_df$phylo]
obs_df$Elevation <- Elevation[obs_df$Site]
obs_df$Temperature <- Temperature[obs_df$Site]

# Standardize predictors for clean beta estimates
obs_df$Elevation_s <- scale(obs_df$Elevation)[, 1]
obs_df$Body_Mass_s <- scale(obs_df$Body_Mass)[, 1]

# Abundance drops as Temperature exceeds Thermal Tolerance
beta_mismatch <- -0.5
obs_df$Thermal_Mismatch <- obs_df$Temperature - obs_df$Thermal_Tol
log_lambda <- 1.0 +
    beta_mismatch * obs_df$Thermal_Mismatch +
    rnorm(N_obs, 0, 0.5)
obs_df$Abundance <- rpois(N_obs, exp(log_lambda))


# ==========================================
# 4. Model Fitting
# ==========================================
cat("Fitting Dual-Scale because model...\n")

eqs <- list(
    Thermal_Tol ~ Body_Mass_s,
    Temperature ~ Elevation_s + (1 | Site),
    Abundance ~ I(Temperature - Thermal_Tol)
)

fit <- because(
    eqs,
    data = obs_df,
    id_col = "phylo",
    family = c(Abundance = "poisson"),
    structure = tree,
    n.iter = 5000,
    n.burnin = 1000,
    n.cores = 1,
    quiet = TRUE
)

cat("Fitting Naive because model (Ignoring Nested Trees)...\n")

eqs_naive <- list(
    Thermal_Tol ~ Body_Mass_s,
    Temperature ~ Elevation_s,
    Abundance ~ I(Temperature - Thermal_Tol)
)

fit_naive <- because(
    eqs_naive,
    data = obs_df,
    id_col = "phylo",
    family = c(Abundance = "poisson"),
    structure = NULL, # Dropping the phylogenetic tree!
    n.iter = 5000,
    n.burnin = 1000,
    n.cores = 1,
    quiet = TRUE
)


# ==========================================
# 5. Extract Results
# ==========================================
samples <- as.matrix(fit$samples)
samples_naive <- as.matrix(fit_naive$samples)

cat("\n\n--- RESULTS FOR PAPER ---\n")
cat("=== DUAL SCALE (CORRECT) ===\n")
cat(sprintf(
    "Beta Thermal_Tol ~ Body_Mass (True: 0.8): %.2f +/- %.2f\n",
    mean(samples[, "beta_Thermal_Tol_Body_Mass_s"]),
    sd(samples[, "beta_Thermal_Tol_Body_Mass_s"])
))
cat(sprintf(
    "Beta Temperature ~ Elevation (True: -1.5): %.2f +/- %.2f\n",
    mean(samples[, "beta_Temperature_Elevation_s"]),
    sd(samples[, "beta_Temperature_Elevation_s"])
))
cat(sprintf(
    "Beta Abundance ~ Mismatch (True: -0.5): %.2f +/- %.2f\n",
    mean(samples[, "beta_Abundance_Temperature_minus_Thermal_Tol"]),
    sd(samples[, "beta_Abundance_Temperature_minus_Thermal_Tol"])
))

cat("\n=== NAIVE MODEL (IGNORING SPACE/PHYLOGENY) ===\n")
cat(sprintf(
    "Beta Thermal_Tol ~ Body_Mass (True: 0.8): %.2f +/- %.2f\n",
    mean(samples_naive[, "beta_Thermal_Tol_Body_Mass_s"]),
    sd(samples_naive[, "beta_Thermal_Tol_Body_Mass_s"])
))
cat(sprintf(
    "Beta Temperature ~ Elevation (True: -1.5): %.2f +/- %.2f\n",
    mean(samples_naive[, "beta_Temperature_Elevation_s"]),
    sd(samples_naive[, "beta_Temperature_Elevation_s"])
))
cat(sprintf(
    "Beta Abundance ~ Mismatch (True: -0.5): %.2f +/- %.2f\n",
    mean(samples_naive[, "beta_Abundance_Temperature_minus_Thermal_Tol"]),
    sd(samples_naive[, "beta_Abundance_Temperature_minus_Thermal_Tol"])
))

saveRDS(fit, "simulated_dual_scale.rds")
cat("\nSimulation and fitting complete.\n")
