library(because)
library(because.phybase)
library(ape)
library(lme4)

set.seed(42)

# --- 1. Simulate Full Data (Same as 4-Level Community) ---
N_Species <- 50
species_tree <- rtree(N_Species)
species_tree$tip.label <- paste0("Sp_", 1:N_Species)
# Normalize tree depth
species_tree$edge.length <- species_tree$edge.length /
    max(node.depth.edgelength(species_tree))

# Species traits
Body_Mass_s <- rnorm(N_Species)
Metabolic_Rate <- 0.8 * Body_Mass_s + rnorm(N_Species, 0, 0.1)
Thermal_Tol <- 0.8 * Metabolic_Rate + rnorm(N_Species, 0, 0.1)
d_species <- data.frame(
    Species = species_tree$tip.label,
    Body_Mass_s,
    Metabolic_Rate,
    Thermal_Tol
)

# Site traits
N_Site <- 15
Elevation_s <- rnorm(N_Site)
U_Resource <- 0.8 * Elevation_s + rnorm(N_Site, 0, 0.5)
Temperature <- -0.6 * Elevation_s + rnorm(N_Site, 0, 0.3)
site_re <- rnorm(N_Site, 0, 0.3)
d_site <- data.frame(Site = paste0("Site_", 1:N_Site), Elevation_s, Temperature)
site_coords <- cbind(runif(N_Site, 0, 10), runif(N_Site, 0, 10))
rownames(site_coords) <- d_site$Site
dist_mat <- as.matrix(dist(site_coords))
site_distances <- exp(-dist_mat / 0.5)

# Survey traits
N_Surveys_per_Site <- 2
N_Surveys <- N_Site * N_Surveys_per_Site
d_survey <- data.frame(
    Survey = paste0("Survey_", 1:N_Surveys),
    Site = rep(d_site$Site, each = N_Surveys_per_Site),
    Wind_Speed = rnorm(N_Surveys, 0, 1)
)
survey_re <- rnorm(N_Surveys, 0, 0.2)
d_survey$survey_total_re <- site_re[match(d_survey$Site, d_site$Site)] +
    survey_re

# Observations
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

# Generate Abundance with true beta = -0.4 for Thermal Mismatch
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

# --- 2. Induce Missing Data ---
# Missing completely at random (MCAR): 20% of species are missing Thermal_Tolerance
missing_species_idx <- sample(1:N_Species, size = floor(0.2 * N_Species))
missing_species <- d_species$Species[missing_species_idx]

d_species_miss <- d_species
d_species_miss$Thermal_Tol[missing_species_idx] <- NA

d_obs_miss <- d_obs
d_obs_miss$Thermal_Tol[d_obs_miss$Species %in% missing_species] <- NA

cat(sprintf(
    "Induced missingness: %d out of %d observations are missing Thermal_Tol\n",
    sum(is.na(d_obs_miss$Thermal_Tol)),
    nrow(d_obs_miss)
))

# --- 3. Fit PGLS / GLMM (Listwise Deletion) ---
# GLMM requires dropping the NA rows. It can't impute the missing traits.
d_obs_complete <- na.omit(d_obs_miss)

cat(sprintf(
    "\nFitting GLMM on complete cases (%d rows retained)...\n",
    nrow(d_obs_complete)
))
# Center Wind_Speed for GLMM stability
d_obs_complete$Wind_Speed_c <- scale(d_obs_complete$Wind_Speed)
d_obs_complete$Thermal_Mismatch <- d_obs_complete$Temperature -
    d_obs_complete$Thermal_Tol

# glmer model roughly approximating the structural equations (lacking phylogeny/spatial)
glmm_fit <- glmer(
    Abundance ~ scale(Foraging_Intensity) +
        Wind_Speed_c +
        Thermal_Mismatch +
        (1 | Site) +
        (1 | Survey) +
        (1 | Species),
    data = d_obs_complete,
    family = poisson
)
glmm_summary <- summary(glmm_fit)
print(glmm_summary$coefficients[
    "Thermal_Mismatch",
    c("Estimate", "Std. Error")
])

# --- 4. Fit because (Full Information Joint Imputation) ---
data_list_miss <- list(
    species = d_species_miss,
    site = d_site,
    survey = d_survey[, c("Survey", "Site", "Wind_Speed")],
    obs = d_obs_miss[, c(
        "Survey",
        "Species",
        "Foraging_Intensity",
        "Abundance"
    )]
)

structures <- list(
    phylo = species_tree,
    spatial = site_distances
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

cat("\nFitting Joint Bayesian Model (because) with Imputation...\n")
bsem_fit <- because(
    eqs,
    data = data_list_miss,
    levels = list(
        species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
        site = c("Elevation_s", "Temperature", "U_Resource", "Site"),
        survey = c("Wind_Speed", "Survey"),
        obs = c("Foraging_Intensity", "Abundance")
    ),
    hierarchy = "site > survey > obs; species > obs",
    link_vars = list(site = "Site", survey = "Survey", species = "Species"),
    latent = "U_Resource",
    family = c(Abundance = "poisson"),
    structure = structures,
    n.iter = 5000,
    n.burnin = 2000,
    n.cores = 1,
    quiet = FALSE
)

cat("\n--- True Beta for Mismatch: -0.40 ---\n")
cat("GLMM (Listwise Deletion) Mismatch Beta:\n")
print(glmm_summary$coefficients["Thermal_Mismatch", ])

cat("\nBecause (Imputation) Mismatch Beta:\n")
# The summary from because() is a list with 'statistics' and 'quantiles'
sum_stats <- bsem_fit$summary$statistics
mismatch_row <- grep("Temperature_minus_Thermal_Tol", rownames(sum_stats))
if (length(mismatch_row) > 0) {
    print(sum_stats[mismatch_row, ])
} else {
    print("Parameter 'Temperature_minus_Thermal_Tol' not found in summary.")
}
