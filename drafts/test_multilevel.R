devtools::load_all(".")
devtools::load_all("../because.phybase")
library(ape)
library(MASS)
set.seed(42)

N_species <- 50
tree <- rcoal(N_species)
tree$tip.label <- paste0("Sp_", 1:N_species)
vcv_mat <- vcv(tree)

lambda_true <- 0.7
phylo_effects <- mvrnorm(
    1,
    mu = rep(0, N_species),
    Sigma = lambda_true * vcv_mat
)
names(phylo_effects) <- tree$tip.label
Body_Mass <- rnorm(N_species, 10, 2)
names(Body_Mass) <- tree$tip.label

# Species level trait
Thermal_Tol <- 0.5 * Body_Mass + phylo_effects + rnorm(N_species, 0, 0.1)
names(Thermal_Tol) <- tree$tip.label

N_sites <- 50
Elevation <- rnorm(N_sites, 1000, 300)
names(Elevation) <- paste0("Site_", 1:N_sites)
Temperature <- -0.01 * Elevation + rnorm(N_sites, 0, 2)
names(Temperature) <- paste0("Site_", 1:N_sites)

N_obs <- 500
obs_df <- data.frame(
    phylo = sample(tree$tip.label, N_obs, replace = TRUE),
    Site = sample(paste0("Site_", 1:N_sites), N_obs, replace = TRUE)
)

# Populate duplicated data
obs_df$Body_Mass <- Body_Mass[obs_df$phylo]
obs_df$Thermal_Tol <- Thermal_Tol[obs_df$phylo]
obs_df$Elevation <- Elevation[obs_df$Site]
obs_df$Temperature <- Temperature[obs_df$Site]

# Observation level trait
# Abundance drops if Temperature > Thermal_Tol -> Thermal Mismatch
obs_df$Thermal_Mismatch <- obs_df$Temperature - obs_df$Thermal_Tol
# We'll just define log_lambda as depending on Thermal_Mismatch
log_lambda <- 2.0 - 0.2 * obs_df$Thermal_Mismatch + rnorm(N_obs, 0, 0.5)
obs_df$Abundance <- rpois(N_obs, exp(log_lambda))

eqs <- list(
    Thermal_Tol ~ Body_Mass,
    Temperature ~ Elevation,
    Abundance ~ I(Temperature - Thermal_Tol)
)

fit <- because(
    eqs,
    data = obs_df,
    id_col = "phylo",
    family = c(Abundance = "poisson"),
    structure = tree,
    n.iter = 1000,
    n.burnin = 500,
    quiet = FALSE
)

print(fit$summary)
