
# ============================================================
# SERVER BENCHMARK: because vs brms vs phylosem vs MCMCglmm
# ============================================================
# Synchronized with Section 5 Case Study (Phil. Trans. B)
# Includes: Hierarchy + REAL Phylogeny + REAL Spatial Autocorrelation

# 1. SETUP & INSTALLATION
options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("devtools", "dplyr", "ape", "brms", "MASS", "Rcpp", 
              "nimble", "phylolm", "phylopath", "MCMCglmm")
lapply(packages, function(p) if(!require(p, character.only=TRUE)) install.packages(p))

cat("\n--- Installing 'because' Ecosystem ---\n")
devtools::install_github("achazhardenberg/because", upgrade = "never", force = TRUE)
devtools::install_github("achazhardenberg/because.phybase", upgrade = "never", force = TRUE)

library(because)
library(because.phybase)
library(brms)
library(ape)
library(dplyr)
library(nimble)
library(MCMCglmm)

n_cores <- parallel::detectCores()

# 2. DATA GENERATION (Hierarchical + REAL Phylo + REAL Spatial)
set.seed(42)
N_Species <- 50
N_Site <- 30
N_Surveys_per_Site <- 3

# --- Space & Geography ---
site_coords <- cbind(runif(N_Site, 0, 10), runif(N_Site, 0, 10))
rownames(site_coords) <- paste0("Site_", 1:N_Site)
V_spatial <- exp(-as.matrix(dist(site_coords)) / 1.0) 

# --- REAL Phylogeny ---
species_tree <- compute.brlen(rtree(N_Species))
species_tree$tip.label <- paste0("Sp_", 1:N_Species)
Vcv_phylo <- vcv.phylo(species_tree)
# Pagel's lambda for all trait residuals
V_lambda  <- 0.60 * (Vcv_phylo/max(Vcv_phylo)) + 0.40 * diag(N_Species)

# --- REAL Spatial Confounding ---
Elevation_s <- as.vector(MASS::mvrnorm(1, rep(0, N_Site), V_spatial))

# --- Local Environment ---
U_Resource   <- 0.80 * Elevation_s + rnorm(N_Site, 0, 0.40)
NDVI         <- 0.50 * U_Resource + 0.30 * Elevation_s + rnorm(N_Site, 0, 0.1)
Flower_Cover <- 0.60 * U_Resource + 0.20 * Elevation_s + rnorm(N_Site, 0, 0.1)
d_site <- data.frame(Site = rownames(site_coords), Elevation_s, NDVI, Flower_Cover)

# --- REAL Evolutionary Confounding (Weighted Residuals) ---
Body_Mass_s    <- as.vector(MASS::mvrnorm(1, rep(0, N_Species), V_lambda))
eps_MR          <- as.vector(MASS::mvrnorm(1, rep(0, N_Species), 0.15^2 * V_lambda))
Metabolic_Rate <- 0.80 * Body_Mass_s + eps_MR
eps_TT          <- as.vector(MASS::mvrnorm(1, rep(0, N_Species), 0.10^2 * V_lambda))
Thermal_Tol    <- 0.60 * Metabolic_Rate + 0.20 * Body_Mass_s + eps_TT
d_species      <- data.frame(Species = species_tree$tip.label, Body_Mass_s, Metabolic_Rate, Thermal_Tol)

# --- Observations ---
site_of_survey  <- rep(1:N_Site, each = N_Surveys_per_Site)
Temperature     <- -1.50 * Elevation_s[site_of_survey] + rnorm(N_Site * N_Surveys_per_Site, 0, 0.2)
Wind_Speed      <- rnorm(N_Site * N_Surveys_per_Site, 0, 1)
d_survey        <- data.frame(Survey = paste0("Survey_", 1:length(Temperature)), 
                              Site = d_site$Site[site_of_survey], Temperature, Wind_Speed)

d_obs <- expand.grid(Survey = d_survey$Survey, Species = d_species$Species, stringsAsFactors = FALSE) %>%
    inner_join(d_survey, by = "Survey") %>%
    inner_join(d_species, by = "Species") %>%
    inner_join(d_site, by = "Site") %>%
    mutate(
        log_lambda = 0.50 * Flower_Cover - 0.50 * (Temperature - Thermal_Tol) + rnorm(n(), 0, 0.10),
        Abundance = rpois(n(), exp(log_lambda))
    ) %>%
    mutate_if(is.character, as.factor)

# 3. FIT MODELS
structures <- list(phylo = species_tree, spatial = V_spatial)
eqs <- list(NDVI ~ U_Resource + Elevation_s, Flower_Cover ~ U_Resource + Elevation_s,
            Metabolic_Rate ~ Body_Mass_s, Thermal_Tol ~ Metabolic_Rate + Body_Mass_s,
            Abundance ~ Flower_Cover + NDVI + Body_Mass_s + Wind_Speed + I(Temperature - Thermal_Tol) + (1|Site) + (1|Survey))

# --- because (The Simultaneous Master - JAGS) ---
cat("\n--- Running because (JAGS with REAL Spatial & REAL Phylo) ---\n")
fit_because_jags <- because(eqs, data = d_obs, structure = structures,
                       levels = list(species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
                                     site    = c("Elevation_s", "NDVI", "Flower_Cover", "Site"),
                                     survey  = c("Temperature", "Wind_Speed", "Survey"),
                                     obs     = c("Abundance")),
                       hierarchy = "site > survey > obs; species > obs",
                       link_vars = list(site = "Site", survey = "Survey", species = "Species"),
                       latent = "U_Resource", latent_method = "explicit", 
                       engine = "jags", parallel = TRUE)
saveRDS(fit_because_jags, "server_because_jags_final_fit.rds")

# --- because (The Simultaneous Master - NIMBLE) ---
cat("\n--- Running because (NIMBLE with REAL Spatial & REAL Phylo) ---\n")
fit_because_nimble <- because(eqs, data = d_obs, structure = structures,
                       levels = list(species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
                                     site    = c("Elevation_s", "NDVI", "Flower_Cover", "Site"),
                                     survey  = c("Temperature", "Wind_Speed", "Survey"),
                                     obs     = c("Abundance")),
                       hierarchy = "site > survey > obs; species > obs",
                       link_vars = list(site = "Site", survey = "Survey", species = "Species"),
                       latent = "U_Resource", latent_method = "explicit", 
                       engine = "nimble", parallel = TRUE)
saveRDS(fit_because_nimble, "server_because_nimble_final_fit.rds")

# --- brms (The Expert Bayesian Baseline) ---
cat("\n--- Running brms (Spatial + Phylo PGLMM) ---\n")
# Fair Comparison: Give brms all species-level traits
bf_abund <- bf(Abundance ~ Flower_Cover + NDVI + Body_Mass_s + Metabolic_Rate + Thermal_Tol + (1 | gr(Site, cov = V_spatial)) + (1 | gr(Species, cov = Vcv_phylo)))

fit_brms <- brm(bf_abund, data = d_obs, data2 = list(V_spatial = V_spatial, Vcv_phylo = Vcv_phylo), 
                family = poisson(), cores = n_cores, iter = 4000)
saveRDS(fit_brms, "server_brms_final_fit.rds")

# --- MCMCglmm (The Sequential Spatial PGLMM) ---
cat("\n--- Running MCMCglmm (Spatial PGLMM) ---\n")
inv_phylo   <- inverseA(species_tree, nodes = "tips", scale = TRUE)$Ainv
inv_spatial <- as(solve(V_spatial), "dgCMatrix")
rownames(inv_spatial) <- colnames(inv_spatial) <- rownames(V_spatial)

prior_mcmc <- list(R = list(V = 1, nu = 0.002), 
                   G = list(G1 = list(V = 1, nu = 0.002), 
                            G2 = list(V = 1, nu = 0.002),
                            G3 = list(V = 1, nu = 0.002)))

fit_mcmc <- MCMCglmm(Abundance ~ Flower_Cover + NDVI + Body_Mass_s + Wind_Speed + Temperature + Thermal_Tol,
                     random = ~Site + Survey + Species, 
                     ginverse = list(Species = inv_phylo, Site = inv_spatial), 
                     prior = prior_mcmc, family = "poisson", data = d_obs,
                     nitt = 60000, burnin = 10000)
saveRDS(fit_mcmc, "server_mcmcglmm_final_fit.rds")

# FINAL SAVE
save.image("benchmark_FINAL_results.RData")
cat("\n[Success] Benchmark complete.\n")
