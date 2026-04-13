
# ============================================================
# SERVER BENCHMARK: because vs brms vs phylosem
# ============================================================
# This script is designed for aservers with multiple cores.
# It installs dependencies, runs the full suite, and saves findings.

# 1. SETUP & INSTALLATION
if (!require("devtools")) install.packages("devtools")
if (!require("dplyr")) install.packages("dplyr")
if (!require("ape")) install.packages("ape")
if (!require("brms")) install.packages("brms")
if (!require("MASS")) install.packages("MASS")
if (!require("Rcpp")) install.packages("Rcpp")

# Install 'because' ecosystem from GitHub
cat("\n--- Installing 'because' Ecosystem from GitHub ---\n")
devtools::install_github("achazhardenberg/because", upgrade = "never", force = TRUE)
devtools::install_github("achazhardenberg/because.phybase", upgrade = "never", force = TRUE)

library(because)
library(because.phybase)
library(brms)
library(ape)
library(dplyr)
library(phylosem) # Assumes phylosem is already available or installed via devtools

# Configuration
n_cores <- parallel::detectCores()
cat("\nDetected Cores:", n_cores, "\n")

# ============================================================
# SECTION A: DATA GENERATION (Identical to Manuscript)
# ============================================================
set.seed(42)
N_Species <- 50
N_Site <- 30
N_Surveys_per_Site <- 3
N_Surveys <- N_Site * N_Surveys_per_Site

species_tree <- rtree(N_Species)
species_tree$tip.label <- paste0("Sp_", 1:N_Species)
Vcv <- vcv.phylo(species_tree)
V_lambda <- 0.60 * Vcv + (1 - 0.60) * diag(diag(Vcv))

Body_Mass_s    <- as.vector(MASS::mvrnorm(1, rep(0, N_Species), V_lambda))
Metabolic_Rate <- 0.80 * Body_Mass_s + as.vector(MASS::mvrnorm(1, rep(0, N_Species), 0.15^2 * V_lambda))
Thermal_Tol    <- 0.60 * Metabolic_Rate + 0.20 * Body_Mass_s + as.vector(MASS::mvrnorm(1, rep(0, N_Species), 0.10^2 * V_lambda))
d_species <- data.frame(Species = species_tree$tip.label, Body_Mass_s, Metabolic_Rate, Thermal_Tol)

Elevation_s  <- rnorm(N_Site)
U_Resource   <- 0.80 * Elevation_s + rnorm(N_Site, 0, 0.50)
NDVI         <- 0.50 * U_Resource + 0.30 * Elevation_s + rnorm(N_Site, 0, 0.20)
Flower_Cover <- 0.60 * U_Resource + 0.20 * Elevation_s + rnorm(N_Site, 0, 0.20)
d_site <- data.frame(Site = paste0("Site_", 1:N_Site), Elevation_s, NDVI, Flower_Cover)

site_of_survey  <- rep(1:N_Site, each = N_Surveys_per_Site)
Temperature     <- -1.50 * Elevation_s[site_of_survey] + rnorm(N_Surveys, 0, 0.30)
Wind_Speed      <- rnorm(N_Surveys, 0, 1)
d_survey        <- data.frame(Survey = paste0("Survey_", 1:N_Surveys), Site = d_site$Site[site_of_survey], Temperature, Wind_Speed)

d_obs <- expand.grid(Survey = d_survey$Survey, Species = d_species$Species, stringsAsFactors = FALSE) %>%
    inner_join(d_survey, by = "Survey") %>%
    inner_join(d_species, by = "Species") %>%
    inner_join(d_site, by = "Site") %>%
    mutate(
        mismatch   = Temperature - Thermal_Tol,
        log_lambda = 0.50 * Flower_Cover - 0.20 * Wind_Speed - 0.50 * mismatch + rnorm(n(), 0, 0.1),
        Abundance   = rpois(n(), exp(log_lambda)),
        is_species = !duplicated(Species),
        is_site    = !duplicated(Site),
        is_survey  = !duplicated(Survey)
    )

# ============================================================
# SECTION B: "because" ANALYSIS
# ============================================================
cat("\n--- Fitting 'because' (Parallel) ---\n")
start_time <- Sys.time()
fit_because <- because(
    Abundance ~ Flower_Cover + Wind_Speed + mismatch,
    Thermal_Tol ~ Metabolic_Rate + Body_Mass_s,
    Metabolic_Rate ~ Body_Mass_s,
    Temperature ~ Elevation_s,
    NDVI ~ U_Resource + Elevation_s,
    Flower_Cover ~ U_Resource + Elevation_s,
    data = d_obs,
    hierarchy = "site > survey > obs; species > obs",
    link_vars = list(site = "Site", survey = "Survey", species = "Species"),
    latent    = "U_Resource",
    engine    = "jags", parallel = TRUE, 
    n.iter    = 12500, n.burnin = 2500
)
time_because <- difftime(Sys.time(), start_time, units="secs")
saveRDS(fit_because, "server_because_fit.rds")

# ============================================================
# SECTION C: "brms" EXPERT (Resolution-Aware)
# ============================================================
cat("\n--- Fitting 'brms' Expert (Parallel) ---\n")
bf_abund <- bf(Abundance | subset(is_survey) ~ mi(Flower_Cover) + Wind_Speed + mi(Temperature) + mi(Thermal_Tol))
bf_flower <- bf(Flower_Cover | subset(is_site) ~ mi(U_Resource) + Elevation_s)
bf_temp  <- bf(Temperature | subset(is_survey) ~ Elevation_s)
bf_tt    <- bf(Thermal_Tol | subset(is_species) ~ mi(Metabolic_Rate) + Body_Mass_s)
bf_mr    <- bf(Metabolic_Rate | subset(is_species) ~ Body_Mass_s)
bf_u     <- bf(U_Resource | subset(is_site) ~ Elevation_s)

fit_brms_expert <- brm(
    bf_abund + bf_flower + bf_temp + bf_tt + bf_mr + bf_u + set_rescor(FALSE),
    data = d_obs, family = gaussian(), 
    cores = n_cores, iter = 1200, warmup = 400,
    control = list(adapt_delta = 0.95)
)
saveRDS(fit_brms_expert, "server_brms_expert_fit.rds")

# ============================================================
# SECTION D: "phylosem"
# ============================================================
cat("\n--- Fitting 'phylosem' ---\n")
d_collapsed <- d_obs %>%
    group_by(Species) %>%
    summarize(
        Abundance_mean = mean(Abundance),
        Flower_mean    = mean(Flower_Cover),
        NDVI_mean      = mean(NDVI),
        Wind_mean      = mean(Wind_Speed),
        mismatch_mean  = mean(mismatch),
        .groups = "drop"
    ) %>%
    inner_join(d_species, by = "Species")

d_collapsed$U_Resource <- NA_real_
model_syntax <- "
    U_Resource =~ 1*Flower_mean + NDVI_mean
    Abundance_mean ~ U_Resource + Wind_mean + mismatch_mean
    Metabolic_Rate ~ Body_Mass_s
    Thermal_Tol ~ Metabolic_Rate + Body_Mass_s
"
fit_psem <- phylosem(sem = model_syntax, data = d_collapsed, 
                     tree = species_tree, family = "poisson")
saveRDS(fit_psem, "server_psem_fit.rds")

# FINAL SAVE
save.image("benchmark_server_results.RData")
cat("\n[Success] Benchmark complete. Results saved to 'benchmark_server_results.RData'.\n")
