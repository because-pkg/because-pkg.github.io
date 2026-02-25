devtools::load_all(".")
devtools::load_all("../because.phybase")
library(ape)
library(MASS)

set.seed(42)

# 1. Evolutionary level: Phylogeny (50 species)
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

# 2. Ecological level: Sites (50 sites)
N_sites <- 50
site_coords <- matrix(runif(N_sites * 2, 0, 100), ncol = 2)
rownames(site_coords) <- paste0("Site_", 1:N_sites)
dist_mat <- as.matrix(dist(site_coords))
V_spatial <- exp(-dist_mat / 20) # Spatial correlation matrix

site_elevation <- rnorm(N_sites, mean = 1000, sd = 300)
site_effects_abundance <- mvrnorm(
    1,
    mu = rep(0, N_sites),
    Sigma = 0.5^2 * V_spatial
)
names(site_effects_abundance) <- rownames(site_coords)

# 3. Observations (500 observations)
N_obs <- 500
obs_df <- data.frame(
    phylo = sample(tree$tip.label, N_obs, replace = TRUE),
    Site = sample(rownames(site_coords), N_obs, replace = TRUE)
)

obs_df$Elevation <- site_elevation[as.numeric(gsub("Site_", "", obs_df$Site))]
obs_df$phylo_u <- phylo_effects[obs_df$phylo]
obs_df$site_v <- site_effects_abundance[obs_df$Site]

obs_df$Elevation_s <- scale(obs_df$Elevation)[, 1]

# Simulate Thermal_Tol
beta_1 <- 1.5
obs_df$Thermal_Tol <- beta_1 *
    obs_df$Elevation_s +
    obs_df$phylo_u +
    rnorm(N_obs, 0, 0.5)

# Simulate Abundance
beta_2 <- 0.4
log_lambda <- 0.0 + beta_2 * obs_df$Thermal_Tol + obs_df$site_v
obs_df$Abundance <- rpois(N_obs, exp(log_lambda))

cat("Fitting because model...\n")

eqs <- list(
    Thermal_Tol ~ Elevation_s,
    Abundance ~ Thermal_Tol + (1 | Site)
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
    parallel = FALSE,
    verbose = FALSE,
    quiet = FALSE
)

# Extract True vs Estimated params
# because() summary method returns formatted text, so we calculate means from raw samples
samples <- as.matrix(fit$samples)

cat("\n\n--- RESULTS FOR PAPER ---\n")
cat(sprintf(
    "Beta Thermal_Tol ~ Elevation (True: 1.5): %.2f\n",
    mean(samples[, "beta_Thermal_Tol_Elevation_s"])
))
cat(sprintf(
    "Beta Abundance ~ Thermal_Tol (True: 0.4): %.2f\n",
    mean(samples[, "beta_Abundance_Thermal_Tol"])
))
cat(sprintf(
    "Phylo Signal Lambda (True: 0.7): %.1f\n",
    mean(samples[, "sigma_Thermal_Tol_phylo"]^2) /
        (mean(samples[, "sigma_Thermal_Tol_phylo"]^2) +
            mean(samples[, "sigma_Thermal_Tol_res"]^2))
))

cat("\nD-Separation Tests:\n")
print(fit$dsep_results)

saveRDS(fit, "simulated_dual_scale.rds")
cat("Simulation and fitting complete.\n")
