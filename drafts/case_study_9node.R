# ============================================================
# Supplementary Code S1
# Case Study: Section 5 of
#   "Because we can: Joint Causal Inference across Evolutionary
#    and Ecological Scales"
#   Phil. Trans. Roy. Soc. B
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
#
# Authors: von Hardenberg & Gonzalez-Voyer
# ============================================================

library(devtools)
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.phybase")
library(ape)
library(MASS)

set.seed(42)

# ============================================================
# SECTION A: GENERATIVE PARAMETERS
# ============================================================

# --- Species-level causal paths ---
BETA_BM_MR          <-  0.80   # Body_Mass_s    → Metabolic_Rate
BETA_MR_TT          <-  0.60   # Metabolic_Rate → Thermal_Tol
BETA_BM_TT          <-  0.20   # Body_Mass_s    → Thermal_Tol (direct)

# --- Site-level causal paths ---
# Elevation_s drives both the long-run mean Temperature (via site-level effect
# transmitted into each survey) and the two remote-sensing variables.
BETA_ELEV_URES      <-  0.80   # Elevation_s  → U_Resource (latent)
BETA_URES_NDVI      <-  0.50   # U_Resource   → NDVI
BETA_ELEV_NDVI      <-  0.30   # Elevation_s  → NDVI  (direct)
BETA_URES_FC        <-  0.60   # U_Resource   → Flower_Cover
BETA_ELEV_FC        <-  0.20   # Elevation_s  → Flower_Cover (direct)

# --- Survey-level causal paths ---
# Temperature is measured at each visit. Elevation sets the site mean;
# survey-to-survey weather variation (sd = 0.30) adds realistic noise.
BETA_ELEV_TEMP      <- -1.50   # Elevation_s  → Temperature (site-level predictor)
SD_TEMP_WITHIN      <-  0.30   # Within-site sd of Temperature across visits

# --- Observation-level causal paths ---
BETA_FC_ABUND       <-  0.50   # Flower_Cover          → log(Abundance)
BETA_WIND_ABUND     <- -0.20   # Wind_Speed            → log(Abundance)
BETA_MISMATCH_ABUND <- -0.50   # I(Temp-Thermal_Tol)  → log(Abundance)

# --- Phylogenetic signal ---
PAGEL_LAMBDA        <-  0.60   # True Pagel's λ for species traits

# ============================================================
# SECTION B: SIMULATION DESIGN
# ============================================================
N_Species          <- 50
N_Site             <- 30   # 30 sites
N_Surveys_per_Site <-  3   # 3 repeat visits per site per season
N_Surveys          <- N_Site * N_Surveys_per_Site   # 90 surveys
# N_obs = 90 surveys × 50 species = 4500

# ============================================================
# SECTION C: SPECIES LEVEL  (S = 50)
# ============================================================
species_tree <- rtree(N_Species)
species_tree$tip.label  <- paste0("Sp_", 1:N_Species)
species_tree$edge.length <- species_tree$edge.length /
    max(node.depth.edgelength(species_tree))

# Pagel's lambda VCV: V_lambda = lambda*Vcv + (1-lambda)*diag(Vcv)
Vcv      <- vcv.phylo(species_tree)
V_lambda <- PAGEL_LAMBDA * Vcv + (1 - PAGEL_LAMBDA) * diag(diag(Vcv))

Body_Mass_s <- as.vector(mvrnorm(1, rep(0, N_Species), V_lambda))
# REAL Evolutionary Confounding: residuals are now weighted
eps_MR      <- as.vector(mvrnorm(1, rep(0, N_Species), 0.15^2 * V_lambda))
Metabolic_Rate <- BETA_BM_MR * Body_Mass_s + eps_MR
eps_TT      <- as.vector(mvrnorm(1, rep(0, N_Species), 0.10^2 * V_lambda))
Thermal_Tol <- BETA_MR_TT * Metabolic_Rate + BETA_BM_TT * Body_Mass_s + eps_TT

names(Body_Mass_s) <- names(Metabolic_Rate) <- names(Thermal_Tol) <-
    species_tree$tip.label

d_species <- data.frame(
    Species        = species_tree$tip.label,
    Body_Mass_s    = Body_Mass_s,
    Metabolic_Rate = Metabolic_Rate,
    Thermal_Tol    = Thermal_Tol
)

cat(sprintf("\nSpecies traits (N = %d, Pagel lambda = %.2f)\n",
            N_Species, PAGEL_LAMBDA))
cat(sprintf("  Body_Mass_s    : mean = %.2f, sd = %.2f\n",
            mean(Body_Mass_s), sd(Body_Mass_s)))
cat(sprintf("  Metabolic_Rate : mean = %.2f, sd = %.2f\n",
            mean(Metabolic_Rate), sd(Metabolic_Rate)))
cat(sprintf("  Thermal_Tol    : mean = %.2f, sd = %.2f\n",
            mean(Thermal_Tol), sd(Thermal_Tol)))

# ============================================================
# SECTION D: SITE LEVEL  (M = 30 sites)
# Fixed attributes — measured once per site across the season.
# ============================================================
# --- REAL Spatial Confounding: Elevation is now Spatially Weighted ---
site_coords <- cbind(runif(N_Site, 0, 10), runif(N_Site, 0, 10))
rownames(site_coords) <- paste0("Site_", 1:N_Site)
V_spatial <- exp(-as.matrix(dist(site_coords)) / 1.0) 

Elevation_s <- as.vector(MASS::mvrnorm(1, rep(0, N_Site), V_spatial))
U_Resource   <- BETA_ELEV_URES * Elevation_s + rnorm(N_Site, 0, 0.50)  # LATENT

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
    Flower_Cover = Flower_Cover,
    U_Resource   = NA_real_   # Latent: JAGS will estimate from NDVI equation
)

# ============================================================
# SECTION E: SURVEY LEVEL  (3 surveys per site × 30 sites = 90 surveys)
# Variables measured at each field visit.
# ============================================================
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

# ============================================================
# SECTION F: OBSERVATION LEVEL  (50 species × 90 surveys = 4500)
# One Pollard walk abundance count per species per survey.
# ============================================================
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

cat(sprintf("\nDataset: %d obs | %d species | %d sites | %d surveys\n",
            nrow(d_obs), N_Species, N_Site, N_Surveys))
cat(sprintf("Mean abundance: %.2f | Zero count: %.1f%%\n",
            mean(d_obs$Abundance), 100 * mean(d_obs$Abundance == 0)))
cat(sprintf("Temperature: mean = %.2f, sd = %.2f (range %.1f to %.1f)\n",
            mean(Temperature), sd(Temperature), min(Temperature), max(Temperature)))

# ============================================================
# SECTION G: PACKAGE DATA STRUCTURES
# ============================================================
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
# SECTION H: STRUCTURAL EQUATIONS  (9-Node MAG)
#
#  Node            Level    Type
#  Elevation_s     Site     Exogenous
#  Body_Mass_s     Species  Exogenous
#  NDVI            Site     ~ U_Resource + Elevation_s
#  Flower_Cover    Site     ~ U_Resource + Elevation_s
#  U_Resource      Site     LATENT — induces NDVI ↔ Flower_Cover in MAG
#  Metabolic_Rate  Species  ~ Body_Mass_s
#  Thermal_Tol     Species  ~ Metabolic_Rate + Body_Mass_s
#  Temperature     Survey   ~ Elevation_s  (site predictor → survey response)
#  Wind_Speed      Survey   Exogenous at survey level
#  Abundance       Obs      ~ Flower_Cover + Wind_Speed + I(Temp-Thermal_Tol)
#
#  Note: only 9 distinct variable nodes; Elevation_s and Body_Mass_s are
#  exogenous (no equation), counted in the DAG but not estimated.
# ============================================================
eqs <- list(
    # --- Site level ---
    NDVI         ~ U_Resource + Elevation_s,
    Flower_Cover ~ U_Resource + Elevation_s,
    # --- Species level ---
    Body_Mass_s    ~ 1,
    Metabolic_Rate ~ Body_Mass_s,
    Thermal_Tol    ~ Metabolic_Rate + Body_Mass_s,
    # --- Survey level: Elevation_s (site) predicts Temperature (survey) ---
    Temperature  ~ Elevation_s,
    # --- Observation level ---
    Abundance ~ Flower_Cover + Wind_Speed +
                I(Temperature - Thermal_Tol) +
                (1 | Site) + (1 | Survey)
)

# ============================================================
# SECTION I: FIT MODELS (Two-Step Causal Pipeline)
# ============================================================

# --- STEP 1: CAUSAL VALIDATION (Basis Set Tests) ---
# We fit the d-separation sub-models to test the MAG topology.
cat("\n[STEP 1] Running Causal Validation (d-separation tests)...\n")
fit_val <- because(
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
    dsep         = TRUE,     # Returns the Basis Set tests
    parallel     = TRUE
)

# --- STEP 2: STRUCTURAL ESTIMATION (Joint Model) ---
# Once validated, we fit the full DAG to retrieve path coefficients.
cat("\n[STEP 2] Running Structural Estimation (Joint MAG)...\n")
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
    parallel     = TRUE
)

# (Optional) NIMBLE Engine for high-speed estimation
# fit_nimble <- because(eqs, data = data_list, ..., engine = "nimble")

# Set 'fit' to fit_jags for results below, but use fit_val for the summary
plot_dag(fit_jags)

# ============================================================
# SECTION J: RESULTS
# ============================================================
cat("\n\n=== M-SEPARATION BASIS TESTS (from STEP 1) ===\n")
s_val <- summary(fit_val)
if (!is.null(s_val$results)) {
    res <- s_val$results
    consistent <- sum(res$LowerCI < 0 & res$UpperCI > 0, na.rm = TRUE)
    cat(sprintf("Tests in basis set:                      25\n"))
    cat(sprintf("Tests executed (cross-hierarchy skipped): %d\n", nrow(res)))
    cat(sprintf("Consistent with independence (CI ∋ 0):   %d\n", consistent))
    cat(sprintf("Rejecting independence (CI ∌ 0):         %d\n",
                nrow(res) - consistent))
    cat("\nDetailed results:\n")
    print(res[, intersect(c("Test", "Parameter", "Estimate", "LowerCI", "UpperCI"),
                          colnames(res))], row.names = FALSE)
    # Visual Causal Validation Plot
    plot_dsep(fit_val)
}

cat("\n=== PARAMETER RECOVERY (from STEP 2) ===\n")
if (!is.null(fit_jags$samples)) {
    mat   <- as.matrix(fit_jags$samples)
    betas <- grep("^beta_", colnames(mat), value = TRUE)
    est   <- apply(mat[, betas, drop = FALSE], 2, mean)
    lo    <- apply(mat[, betas, drop = FALSE], 2, quantile, 0.025)
    hi    <- apply(mat[, betas, drop = FALSE], 2, quantile, 0.975)
    for (b in betas)
        cat(sprintf("  %-52s | Est: %+.2f (95%% BCI: %+.2f, %+.2f)\n",
                    b, est[b], lo[b], hi[b]))
}

# ============================================================
# SECTION K: DIAGNOSTICS & PPC
# ============================================================
cat("\n\n=== GENERATING DIAGNOSTIC PLOTS ===\n")
# 1. Structural Comparison
plot_coef(fit_jags)

# 2. Posterior Predictive Checks (The Causal Chain)
cat("\nRunning PPC for the entire causal chain...\n")
# We check both the final outcome (Abundance) and the intermediate causal drivers
diagnostic_nodes <- c("Metabolic_Rate", "Thermal_Tol", "Flower_Cover", "Temperature", "Abundance")

for (node in diagnostic_nodes) {
    cat(sprintf("  PPC for %s...\n", node))
    # Note: in a real session, these would pop up as separate plot windows or facets
    print(pp_check(fit_jags, resp = node, type = "dens_overlay"))
}

cat("\n=== DONE ===\n")
