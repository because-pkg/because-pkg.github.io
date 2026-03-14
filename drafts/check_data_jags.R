library(because)
library(because.phybase)
library(ape)

set.seed(123)
N_Species <- 50
Body_Mass_s <- rnorm(N_Species)
Metabolic_Rate <- 0.5 * Body_Mass_s + rnorm(N_Species, 0, 0.2)
Thermal_Tol <- 0.6 *
    Metabolic_Rate +
    0.3 * Body_Mass_s +
    rnorm(N_Species, 0, 0.2)
d_species <- data.frame(
    Species = paste0("Sp_", 1:N_Species),
    Body_Mass_s,
    Metabolic_Rate,
    Thermal_Tol
)

N_Site <- 15
Elevation_s <- rnorm(N_Site)
U_Resource <- 0.8 * Elevation_s + rnorm(N_Site, 0, 0.5)
Temperature <- -0.6 * Elevation_s + rnorm(N_Site, 0, 0.3)
site_re <- rnorm(N_Site, 0, 0.3)
d_site <- data.frame(Site = paste0("Site_", 1:N_Site), Elevation_s, Temperature)

N_Surveys_per_Site <- 2
N_Surveys <- N_Site * N_Surveys_per_Site
d_survey <- data.frame(
    Survey = paste0("Survey_", 1:N_Surveys),
    Site = rep(d_site$Site, each = N_Surveys_per_Site),
    Wind_Speed = rnorm(N_Surveys, 10, 2)
)
survey_re <- rnorm(N_Surveys, 0, 0.2)
d_survey$survey_total_re <- site_re[match(d_survey$Site, d_site$Site)] +
    survey_re

d_obs <- expand.grid(Survey = d_survey$Survey, Species = d_species$Species)
d_obs <- merge(
    d_obs,
    d_survey[, c("Survey", "Site", "Wind_Speed", "survey_total_re")],
    by = "Survey"
)
d_obs <- merge(d_obs, d_species, by = "Species")
d_obs <- merge(d_obs, d_site[, c("Site", "Temperature")], by = "Site")
obs_U_Res <- U_Resource[match(d_obs$Site, d_site$Site)]
d_obs$Foraging_Intensity <- 0.5 *
    obs_U_Res -
    0.3 * d_obs$Metabolic_Rate +
    d_obs$survey_total_re +
    rnorm(nrow(d_obs), 0, 0.2)
Abundance_lambda <- exp(
    0.5 *
        obs_U_Res +
        0.3 * d_obs$Foraging_Intensity -
        0.2 * d_obs$Wind_Speed -
        0.4 * (d_obs$Temperature - d_obs$Thermal_Tol) +
        d_obs$survey_total_re +
        rnorm(nrow(d_obs), 0, 0.1)
)
d_obs$Abundance <- rpois(nrow(d_obs), Abundance_lambda)

data_list <- list(
    species = d_species,
    site = d_site,
    survey = d_survey[, c("Survey", "Site", "Wind_Speed")],
    obs = d_obs[, c("Survey", "Species", "Foraging_Intensity", "Abundance")]
)
eqs <- list(
    Temperature ~ Elevation_s,
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

bsem <- because(
    eqs,
    data_list,
    latent = "U_Resource",
    latent_method = "explicit",
    n.iter = 0,
    n.adapt = 0,
    n.burnin = 0,
    levels = list(
        species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
        site = c("Elevation_s", "Temperature", "U_Resource", "Site"),
        survey = c("Wind_Speed", "Survey"),
        obs = c("Foraging_Intensity", "Abundance")
    ),
    hierarchy = "site > survey > obs; species > obs",
    link_vars = list(site = "Site", survey = "Survey", species = "Species")
)

cat("\n----- BSEM Data Lengths -----\n")
if ("data" %in% names(bsem)) {
    cat("Length of group_Site is: ", length(bsem$data$group_Site), "\n")
    cat("Length of group_Survey is: ", length(bsem$data$group_Survey), "\n")
    cat("N_obs is: ", bsem$data$N_obs, "\n")
    cat("N_site is: ", bsem$data$N_site, "\n")
    cat("N_survey is: ", bsem$data$N_survey, "\n")
    cat(
        "group_Site sample: ",
        paste(head(bsem$data$group_Site), collapse = ", "),
        "\n"
    )
    cat(
        "group_Survey sample: ",
        paste(head(bsem$data$group_Survey), collapse = ", "),
        "\n"
    )
}
