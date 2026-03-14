library(because)
library(because.phybase)
library(ape)

set.seed(123)

# --- 1. Simulate Species Traits (N=50) ---
N_Species <- 50
tree <- rtree(N_Species)
tree$tip.label <- paste0("Sp_", 1:N_Species)

# Traits: Body Mass -> Metabolic Rate -> Thermal Tolerance
Body_Mass_s <- rnorm(N_Species)
Metabolic_Rate <- 0.5 * Body_Mass_s + rnorm(N_Species, 0, 0.2)
Thermal_Tol <- 0.6 *
    Metabolic_Rate +
    0.3 * Body_Mass_s +
    rnorm(N_Species, 0, 0.2)

d_species <- data.frame(
    Species = tree$tip.label,
    Body_Mass_s = Body_Mass_s,
    Metabolic_Rate = Metabolic_Rate,
    Thermal_Tol = Thermal_Tol
)

# --- 2. Simulate Site scale (N=15) ---
N_Site <- 15
Elevation_s <- rnorm(N_Site)
U_Resource <- 0.8 * Elevation_s + rnorm(N_Site, 0, 0.5) # LATENT Variable
Temperature <- -0.6 * Elevation_s + rnorm(N_Site, 0, 0.3)
NDVI <- -0.3 * U_Resource + 0.5 * Elevation_s + rnorm(N_Site, 0, 0.2)

# Site random effects (unmeasured site variability)
site_re <- rnorm(N_Site, 0, 0.3)

d_site <- data.frame(
    Site = paste0("Site_", 1:N_Site),
    Elevation_s = Elevation_s,
    Temperature = Temperature,
    NDVI = NDVI
)

# Create a spatial distance matrix for the sites
coords <- matrix(runif(N_Site * 2), ncol = 2)
dist_mat <- as.matrix(dist(coords))
V_spatial <- exp(-dist_mat / 0.5)

# --- 3. Simulate Survey scale (N_surveys = 15 sites * 2 days = 30) ---
N_Surveys_per_Site <- 2
N_Surveys <- N_Site * N_Surveys_per_Site

d_survey <- data.frame(
    Survey = paste0("Survey_", 1:N_Surveys),
    Site = rep(d_site$Site, each = N_Surveys_per_Site),
    Wind_Speed = rnorm(N_Surveys, 0, 1)
)

# Survey random effects (unmeasured variability among days/counts)
survey_re <- rnorm(N_Surveys, 0, 0.2)
# Add site RE to surveys
d_survey$survey_total_re <- site_re[match(d_survey$Site, d_site$Site)] +
    survey_re

# --- 4. Simulate Observation scale (Community counting: All species in all surveys) ---
# N_Obs = 120 surveys * 50 species = 6000
d_obs <- expand.grid(
    Survey = d_survey$Survey,
    Species = d_species$Species
)
d_obs$obs_id <- 1:nrow(d_obs)

# Merge in predictors for simulation
d_obs <- merge(
    d_obs,
    d_survey[, c("Survey", "Site", "Wind_Speed", "survey_total_re")],
    by = "Survey"
)
d_obs <- merge(d_obs, d_species, by = "Species")
d_obs <- merge(d_obs, d_site[, c("Site", "Temperature")], by = "Site")

# Internal latent resource (site-level)
obs_U_Res <- U_Resource[match(d_obs$Site, d_site$Site)]

# Foraging Intensity (depends on Resource and Metabolic Rate + Survey influences)
d_obs$Foraging_Intensity <- 0.5 *
    obs_U_Res -
    0.3 * d_obs$Metabolic_Rate +
    d_obs$survey_total_re +
    rnorm(nrow(d_obs), 0, 0.2)

# Abundance (Poisson)
# depends on Resource, Foraging, Wind, and Thermal Mismatch
Abundance_lambda <- exp(
    0.5 *
        obs_U_Res +
        0.3 * d_obs$Foraging_Intensity -
        0.5 * d_obs$Wind_Speed -
        0.8 * (d_obs$Temperature - d_obs$Thermal_Tol) +
        d_obs$survey_total_re +
        rnorm(nrow(d_obs), 0, 0.1)
)
d_obs$Abundance <- rpois(nrow(d_obs), Abundance_lambda)

# --- 5. Package Data and Fit Model ---
data_list <- list(
    species = d_species,
    site = d_site,
    survey = d_survey[, c("Survey", "Site", "Wind_Speed")],
    obs = d_obs[, c(
        "obs_id",
        "Survey",
        "Species",
        "Foraging_Intensity",
        "Abundance"
    )]
)

structures = list(
    phylo = tree,
    spatial = V_spatial
)

eqs <- list(
    Temperature ~ Elevation_s,
    NDVI ~ U_Resource + Elevation_s,
    Body_Mass_s ~ 1,
    Metabolic_Rate ~ Body_Mass_s,
    Thermal_Tol ~ Metabolic_Rate + Body_Mass_s,
    Foraging_Intensity ~ U_Resource +
        Metabolic_Rate +
        (1 | Site) +
        (1 | Survey),
    Abundance ~ U_Resource +
        Foraging_Intensity +
        Wind_Speed +
        I(Temperature - Thermal_Tol) +
        (1 | Site) +
        (1 | Survey)
)

fit_full <- because(
    eqs,
    data = data_list,
    levels = list(
        species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
        site = c("Elevation_s", "Temperature", "NDVI", "U_Resource", "Site"),
        survey = c("Wind_Speed", "Survey"),
        obs = c("Foraging_Intensity", "Abundance")
    ),
    hierarchy = "site > survey > obs; species > obs",
    link_vars = list(site = "Site", survey = "Survey", species = "Species"),
    latent = "U_Resource",
    latent_method = "explicit",
    family = c(Abundance = "poisson"),
    structure = structures,
    priors = list(
        alpha_Abundance = "dnorm(0, 0.01)",
        beta_Abundance_U_Resource = "dnorm(0, 0.1)",
        beta_Abundance_Foraging_Intensity = "dnorm(0, 0.1)",
        beta_Abundance_Wind_Speed = "dnorm(0, 0.1)",
        beta_Abundance_Temperature_minus_Thermal_Tol = "dnorm(0, 0.1)",
        tau_e_Abundance = "dgamma(3, 3)",
        tau_u_Abundance_Site = "dgamma(2, 2)",
        tau_u_Abundance_Survey = "dgamma(2, 2)"
    ),
    n.iter = 1000,
    n.burnin = 500,
    n.cores = 1,
    quiet = FALSE
)

print(fit_full$summary)
