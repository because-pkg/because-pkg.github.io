# Supplementary Code S1
# Achaz von Hardenberg
# 13 April 2026
#
# 9-Node Hierarchical MAG: Evolutionary × Ecological × Local
#
# Field design:
#   - 30 sites visited 3 times each during one flight season
#   - 50 butterfly species observed at all sites
#   - At each site visit (survey):
#       * Wind_Speed (km/h) and Temperature (°C) recorded by observer
#       * One Pollard walk count per species
#   - At each site (once per season):
#       * NDVI (remote sensing)
#       * Flower_Cover (% nectar plant cover, vegetation subplot)
#       * Elevation_s (GPS/DEM)
#   - Species trait database: Body_Mass_s, Metabolic_Rate, Thermal_Tol
#   - Molecular phylogeny of the community
#
# Totals: 30 sites × 3 surveys = 90 surveys
#         90 surveys × 50 species = 4500 obs


library(because)
library(because.phybase)
library(ape)
library(MASS)

set.seed(42)

# TRUE PARAMETER VALUES FOR SIMULATED DATA

# Species-level causal paths
BETA_BM_MR          <-  0.80   # Body_Mass_s    → Metabolic_Rate
BETA_MR_TT          <-  0.00   # Decoupled
BETA_BM_TT          <-  0.60   # Body_Mass_s    → Thermal_Tol

# Site-level causal paths
# Elevation_s drives the long-run mean Temperature.
# U_Resource (Soil moisture) is an independent latent confounder.
BETA_ELEV_URES      <-  0.00   # Decoupled
BETA_URES_NDVI      <-  0.50   # U_Resource   → NDVI
BETA_ELEV_NDVI      <-  0.00   # Decoupled
BETA_URES_FC        <-  0.60   # U_Resource   → Flower_Cover
BETA_ELEV_FC        <-  0.20   # Elevation_s  → Flower_Cover (direct)

# Survey-level causal paths
# Temperature is measured at each visit. Elevation sets the site mean;
# survey-to-survey weather variation (sd = 0.30) adds realistic noise.
BETA_ELEV_TEMP      <- -1.50   # Elevation_s  → Temperature (site-level predictor)
SD_TEMP_WITHIN      <-  0.30   # Within-site sd of Temperature across visits

#  Observation-level causal paths
BETA_FC_ABUND       <-  0.50   # Flower_Cover          → log(Abundance)
BETA_WIND_ABUND     <- -0.20   # Wind_Speed            → log(Abundance)
BETA_MISMATCH_ABUND <- -0.50   # I(Temp-Thermal_Tol)  → log(Abundance)

# Phylogenetic signal
PAGEL_LAMBDA        <-  0.60   # True Pagel's λ for species traits

# SIMULATION DESIGN

N_Species          <- 50
N_Site             <- 30   # 30 sites
N_Surveys_per_Site <-  3   # 3 repeat visits per site per season
N_Surveys          <- N_Site * N_Surveys_per_Site   # 90 surveys
# N_obs = 90 surveys × 50 species = 4500


# SPECIES LEVEL  (S = 50)

species_tree <- rtree(N_Species)
species_tree$tip.label  <- paste0("Sp_", 1:N_Species)
species_tree$edge.length <- species_tree$edge.length /
    max(node.depth.edgelength(species_tree))

# Pagel's lambda VCV: V_lambda = lambda*Vcv + (1-lambda)*diag(Vcv)
Vcv      <- vcv.phylo(species_tree)
V_lambda <- PAGEL_LAMBDA * Vcv + (1 - PAGEL_LAMBDA) * diag(diag(Vcv))

Body_Mass_s <- as.vector(mvrnorm(1, rep(0, N_Species), V_lambda))
# Evolutionary Confounding: residuals are weighted
eps_MR      <- as.vector(mvrnorm(1, rep(0, N_Species), 0.15^2 * V_lambda))
Metabolic_Rate <- BETA_BM_MR * Body_Mass_s + eps_MR
eps_TT      <- as.vector(mvrnorm(1, rep(0, N_Species), 0.10^2 * V_lambda))
Thermal_Tol <- BETA_BM_TT * Body_Mass_s + eps_TT # No MR effect

names(Body_Mass_s) <- names(Metabolic_Rate) <- names(Thermal_Tol) <-
    species_tree$tip.label

d_species <- data.frame(
    Species        = species_tree$tip.label,
    Body_Mass_s    = Body_Mass_s,
    Metabolic_Rate = Metabolic_Rate,
    Thermal_Tol    = Thermal_Tol
)

# SITE LEVEL  (M = 30 sites)
# Fixed attributes — measured once per site across the season.

# Spatial Confounding: Elevation is  Spatially Weighted
site_coords <- cbind(runif(N_Site, 0, 10), runif(N_Site, 0, 10))
rownames(site_coords) <- paste0("Site_", 1:N_Site)
V_spatial <- exp(-as.matrix(dist(site_coords)) / 1.0)

Elevation_s <- as.vector(MASS::mvrnorm(1, rep(0, N_Site), V_spatial))
# U_Resource is modeled as N(0, 1) standard noise for the MAG demonstration
U_Resource   <- rnorm(N_Site, 0, 1.0)  

# NDVI and Flower_Cover both driven by same latent U_Resource
# → induces the bidirected edge NDVI ↔ Flower_Cover in the MAG
NDVI         <- BETA_URES_NDVI * U_Resource +
                BETA_ELEV_NDVI * Elevation_s +
                rnorm(N_Site, 0, 0.20)
Flower_Cover <- BETA_URES_FC   * U_Resource +
                BETA_ELEV_FC   * Elevation_s +
                rnorm(N_Site, 0, 0.20)

d_site <- data.frame(
    Site         = paste0("Site_", 1:N_Site),
    Elevation_s  = Elevation_s,
    NDVI         = NDVI,
    Flower_Cover = Flower_Cover
)

# SURVEY LEVEL  (3 surveys per site × 30 sites = 90 surveys)
# Variables measured at each field visit.

site_of_survey  <- rep(1:N_Site, each = N_Surveys_per_Site)

# Temperature: site mean set by Elevation + survey-to-survey weather noise
Temperature <- BETA_ELEV_TEMP * Elevation_s[site_of_survey] +
               rnorm(N_Surveys, 0, SD_TEMP_WITHIN)

# Wind_Speed: purely local atmospheric condition at time of visit
Wind_Speed  <- rnorm(N_Surveys, 0, 1)

d_survey <- data.frame(
    Survey      = paste0("Survey_", 1:N_Surveys),
    Site        = d_site$Site[site_of_survey],
    Temperature = Temperature,   # measured per visit with digital thermometer
    Wind_Speed  = Wind_Speed     # measured per visit with anemometer
)

# OBSERVATION LEVEL  (50 species × 90 surveys = 4500)
# One Pollard walk abundance count per survey in which abundance of all 50 species are counted.
d_obs <- expand.grid(
    Survey  = d_survey$Survey,
    Species = d_species$Species,
    stringsAsFactors = FALSE
)
d_obs <- merge(d_obs, d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed")],
               by = "Survey")
d_obs <- merge(d_obs, d_species, by = "Species")
d_obs <- merge(d_obs, d_site[, c("Site", "Flower_Cover")], by = "Site")
d_obs$obs_id <- seq_len(nrow(d_obs))

# Thermal mismatch: temperature at the survey visit vs species thermal optimum
Thermal_Mismatch <- d_obs$Temperature - d_obs$Thermal_Tol

# Abundance: Poisson count
# U_Resource → Abundance ONLY via Flower_Cover (no direct path)
# Metabolic_Rate → Abundance ONLY via Thermal_Tol → I(Temp - Thermal_Tol)
log_lambda <- BETA_FC_ABUND       * d_obs$Flower_Cover   +
              BETA_WIND_ABUND     * d_obs$Wind_Speed      +
              BETA_MISMATCH_ABUND * Thermal_Mismatch      +
              rnorm(nrow(d_obs), 0, 0.10)
d_obs$Abundance <- rpois(nrow(d_obs), exp(log_lambda))



# CREATE LISTS OF DATA STRUCTURES

data_list <- list(
    species = d_species,
    site    = d_site,
    survey  = d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed")],
    obs     = d_obs[, c("obs_id", "Survey", "Species", "Abundance")]
)

structures <- list(
    phylo   = species_tree,   # 50×50 phylogenetic precision (species level)
    spatial = V_spatial        # 30×30 spatial covariance   (site level)
)

# ============================================================

# STRUCTURAL EQUATIONS  (9-Node MAG)

eqs <- list(
    #Site level
    NDVI         ~ U_Resource,
    Flower_Cover ~ U_Resource + Elevation_s,
    # Species level
    Body_Mass_s    ~ 1,
    Metabolic_Rate ~ Body_Mass_s,
    Thermal_Tol    ~ Body_Mass_s,
    # Survey level: Elevation_s (site) predicts Temperature (survey)
    Temperature  ~ Elevation_s,
    # Observation level
    Abundance ~ Flower_Cover + Wind_Speed +
                I(Temperature - Thermal_Tol) +
                (1 | Site) + (1 | Survey)
)

# CONDITIONAL INDEPENDENCIES TESTS
# We fit the d-separation sub-models to test the MAG topology.
# fit_val <- because(
#     eqs,
#     data      = data_list,
#     levels    = list(
#         species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
#         site    = c("Elevation_s", "NDVI", "Flower_Cover", "U_Resource", "Site"),
#         survey  = c("Temperature", "Wind_Speed", "Survey"),
#         obs     = c("Abundance")),
#     hierarchy    = "site > survey > obs; species > obs",
#     link_vars    = list(site = "Site", survey = "Survey", species = "Species"),
#     latent       = "U_Resource",
#     latent_method = "explicit",
#     family       = c(Abundance = "poisson"),
#     structure    = structures,
#     dsep         = TRUE,     # Returns the Basis Set tests
#     parallel     = TRUE,
#     n.iter       = 2500
# )

# FIT FULL JOINT MODEL ---
# Once validated, we fit the full DAG to retrieve path coefficients.
fit_jags <- because(
    eqs,
    data      = data_list,
    levels    = list(
        species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
        site    = c("Elevation_s", "NDVI", "Flower_Cover", "U_Resource", "Site"),
        survey  = c("Temperature", "Wind_Speed", "Survey"),
        obs     = c("Abundance")),
    hierarchy    = "site > survey > obs; species > obs",
    link_vars    = list(site = "Site", survey = "Survey", species = "Species"),
    latent       = "U_Resource",
    latent_method = "explicit",
    family       = c(Abundance = "poisson"),
    structure    = structures,
    dsep         = FALSE,    # Returns the Path Coefficients
    parallel     = TRUE,
    n.cores = 3,
    n.iter = 5000
)

# RESULTS

# M-SEPARATION BASIS TESTS
# summary(fit_val)
# plot_dsep(fit_val)

# PARAMETER RECOVERY
# Inspecting posterior samples for beta coefficients
summary(fit_jags)
plot_coef(fit_jags)

# plot DAG of fitted model
plot_dag(fit_jags)

# POSTERIOR PREDICTIVE CHECKS

# We check both the final outcome (Abundance) and the intermediate causal drivers
diagnostic_nodes <- c("Metabolic_Rate", "Thermal_Tol", "Flower_Cover", "Temperature", "Abundance")

for (node in diagnostic_nodes) {
    # Generating PPC for node...
    print(pp_check(fit_jags, resp = node, type = "dens_overlay"))
}

