set.seed(123)
N_Species <- 50
Body_Mass_s <- rnorm(N_Species)
Metabolic_Rate <- 0.5 * Body_Mass_s + rnorm(N_Species, 0, 0.2)
Thermal_Tol <- 0.6 * Metabolic_Rate + 0.3 * Body_Mass_s + rnorm(N_Species, 0, 0.2)
d_species <- data.frame(Species = paste0("Sp_", 1:N_Species), Body_Mass_s, Metabolic_Rate, Thermal_Tol)

N_Site <- 30
Elevation_s <- rnorm(N_Site)
U_Resource <- 0.8 * Elevation_s + rnorm(N_Site, 0, 0.5)
Temperature <- -0.6 * Elevation_s + rnorm(N_Site, 0, 0.3)
site_re <- rnorm(N_Site, 0, 0.3)
d_site <- data.frame(Site = paste0("Site_", 1:N_Site), Elevation_s, Temperature)

N_Surveys_per_Site <- 1
N_Surveys <- N_Site * N_Surveys_per_Site
d_survey <- data.frame(Survey = paste0("Survey_", 1:N_Surveys), Site = rep(d_site$Site, each = N_Surveys_per_Site), Wind_Speed = rnorm(N_Surveys, 10, 2))
survey_re <- rnorm(N_Surveys, 0, 0.2)
d_survey$survey_total_re <- site_re[match(d_survey$Site, d_site$Site)] + survey_re

d_obs <- expand.grid(Survey = d_survey$Survey, Species = d_species$Species)
d_obs <- merge(d_obs, d_survey[, c("Survey", "Site", "Wind_Speed", "survey_total_re")], by = "Survey")
d_obs <- merge(d_obs, d_species, by = "Species")
d_obs <- merge(d_obs, d_site[, c("Site", "Temperature")], by = "Site")
obs_U_Res <- U_Resource[match(d_obs$Site, d_site$Site)]
d_obs$Foraging_Intensity <- 0.5 * obs_U_Res - 0.3 * d_obs$Metabolic_Rate + d_obs$survey_total_re + rnorm(nrow(d_obs), 0, 0.2)

Abundance_lambda <- exp(
    0.5 * obs_U_Res +
    0.3 * d_obs$Foraging_Intensity -
    0.2 * d_obs$Wind_Speed -
    0.4 * (d_obs$Temperature - d_obs$Thermal_Tol) +
    d_obs$survey_total_re + 
    rnorm(nrow(d_obs), 0, 0.1)
)
print("Lambda summary:")
summary(Abundance_lambda)
