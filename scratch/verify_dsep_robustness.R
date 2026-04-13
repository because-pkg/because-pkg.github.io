# Verification Script v6: Hierarchical D-Sep with Latent
# Achaz von Hardenberg
# 13 April 2026

library(pbapply)
library(parallel)

# Load the local version of because
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")

# 1. Simulate Hierarchical Data (Site / Individual)
set.seed(42)
N_sites <- 5
N_per_site <- 10
N <- N_sites * N_per_site

sites <- rep(1:N_sites, each=N_per_site)
Site_Eff <- rnorm(N_sites)[sites]

# U_Resource is latent (NOT in the site_data)
# But it influences both NDVI and Temperature
U_Resource_vals <- rnorm(N_sites)[sites] 

# Site-level data (NO U_Resource)
site_data <- data.frame(
  Site = 1:N_sites,
  Elevation = rnorm(N_sites)
)

# Observation-level data
obs_data <- data.frame(
  Site = sites,
  NDVI = 0.5 * U_Resource_vals + rnorm(N),
  Temperature = 0.8 * U_Resource_vals + 0.2 * Site_Eff + rnorm(N)
)

# Hierarchical list
data_list <- list(
  site = site_data,
  obs = obs_data
)

# Levels
# U_Resource is in 'levels' for 'site' but NOT in the site_data frame
levels_list <- list(
  site = c("Elevation", "U_Resource", "Site"), 
  obs  = c("NDVI", "Temperature")
)

# Equations (DAG: Elevation -> Temperature; U_Resource -> NDVI, U_Resource -> Temperature)
# Basis set should include Elevation _||_ NDVI | {} (conditional on nothing or Site?)
eqs <- list(
  Temperature ~ Elevation + U_Resource,
  NDVI ~ U_Resource
)

message("\n--- Testing Latent Hierarchical Propagation ---")
message("U_Resource is latent and lives at 'site' level.")

# This should work now without "Variables not found in site dataset: U_Resource"
fit <- because(
  eqs,
  data = data_list,
  levels = levels_list,
  hierarchy = "site > obs",
  link_vars = list(site = "Site"),
  latent = "U_Resource",
  dsep = TRUE,
  parallel = TRUE,
  n.cores = 2,
  n.iter = 500
)

message("\nVerification Complete.")
print(fit$dsep)
