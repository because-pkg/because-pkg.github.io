
library(because)
library(ape)
library(MASS)

set.seed(42)
N_Species <- 10
N_Site <- 10
N_Surveys_per_Site <- 2
N_Surveys <- N_Site * N_Surveys_per_Site

# Minimal data generation
U_Resource <- rnorm(N_Site, 0, 1)
NDVI <- 0.6 * U_Resource + rnorm(N_Site, 0, 0.2)
Flower_Cover <- 0.5 * U_Resource + rnorm(N_Site, 0, 0.2)
Site <- paste0("Site_", 1:N_Site)
d_site <- data.frame(Site = Site, NDVI = NDVI, Flower_Cover = Flower_Cover)

data_list <- list(
    site = d_site
)

eqs <- list(
    NDVI ~ U_Resource,
    Flower_Cover ~ U_Resource
)

fit <- because(
    eqs,
    data = data_list,
    levels = list(site = c("NDVI", "Flower_Cover", "U_Resource", "Site")),
    latent = "U_Resource",
    latent_method = "explicit",
    n.iter = 100,
    n.burnin = 10,
    parallel = FALSE
)

cat("\n--- JAGS MODEL STRING ---\n")
cat(fit$model_string)
cat("\n-------------------------\n")
