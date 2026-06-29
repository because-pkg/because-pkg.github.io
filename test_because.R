library(because)
library(because.phybase)
set.seed(42)
BETA_BM_MR <- 0.80; BETA_BM_TT <- 0.60; BETA_ELEV_FC <- 0.20; BETA_URES_FC <- 0.60; BETA_ELEV_TEMP <- -1.50
BETA_FC_ABUND <- 0.50; BETA_WIND_ABUND <- -0.20; BETA_MISMATCH_ABUND <- -0.50
N_Species <- 50; N_Site <- 30; N_Surveys_per_Site <- 3; N_Surveys <- N_Site * N_Surveys_per_Site
species_tree <- ape::rtree(N_Species)
Vcv <- ape::vcv(species_tree)
Body_Mass_s <- as.vector(MASS::mvrnorm(1, rep(0, N_Species), Vcv))
Thermal_Tol <- BETA_BM_TT * Body_Mass_s + as.vector(MASS::mvrnorm(1, rep(0, N_Species), 0.40^2 * Vcv))
d_species <- data.frame(Species = paste0("Sp_", 1:N_Species), Body_Mass_s = Body_Mass_s, Thermal_Tol = Thermal_Tol)
Elevation_s <- rnorm(N_Site); Flower_Cover <- BETA_ELEV_FC * Elevation_s + rnorm(N_Site)
d_site <- data.frame(Site = paste0("Site_", 1:N_Site), Elevation_s = Elevation_s, Flower_Cover = Flower_Cover)
site_of_survey <- rep(1:N_Site, each = N_Surveys_per_Site)
Temperature <- BETA_ELEV_TEMP * Elevation_s[site_of_survey] + rnorm(N_Surveys)
Wind_Speed <- rnorm(N_Surveys)
d_survey <- data.frame(Survey = paste0("Survey_", 1:N_Surveys), Site = d_site$Site[site_of_survey], Temperature = Temperature, Wind_Speed = Wind_Speed)
d_obs <- expand.grid(Survey = d_survey$Survey, Species = d_species$Species, stringsAsFactors = FALSE)
d_obs <- merge(d_obs, d_survey, by="Survey"); d_obs <- merge(d_obs, d_species, by="Species"); d_obs <- merge(d_obs, d_site, by="Site")
Thermal_Mismatch <- d_obs$Temperature - d_obs$Thermal_Tol
log_lambda <- BETA_FC_ABUND * d_obs$Flower_Cover + BETA_WIND_ABUND * d_obs$Wind_Speed + BETA_MISMATCH_ABUND * Thermal_Mismatch + rnorm(nrow(d_obs), 0, 0.1)
d_obs$Abundance <- rpois(nrow(d_obs), exp(log_lambda))
d_obs$obs_id <- 1:nrow(d_obs)

structures <- list(phylo = species_tree)
data_list <- list(
    species = d_species,
    site    = d_site,
    survey  = d_survey[, c("Survey", "Site", "Temperature", "Wind_Speed")],
    obs     = d_obs[, c("obs_id", "Survey", "Species", "Abundance")]
)
eqs <- list(
    Abundance ~ Body_Mass_s + Thermal_Tol + (1|Site) + (1|Survey) + (1|Species)
)
fit_test <- because(
    eqs,
    data = data_list,
    levels = list(
        species = c("Body_Mass_s", "Thermal_Tol", "Species"),
        site = c("Elevation_s", "Flower_Cover", "Site"),
        survey = c("Temperature", "Wind_Speed", "Survey"),
        obs = c("Abundance")
    ),
    hierarchy = "site > survey > obs; species > obs",
    link_vars = list(site = "Site", survey = "Survey", species = "Species"),
    family = c(Abundance = "negbinomial"),
    structure = structures,
    dsep = FALSE,
    engine = "numpyro"
)
summary(fit_test)
