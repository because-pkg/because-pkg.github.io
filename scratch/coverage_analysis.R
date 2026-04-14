################################################################################
# coverage_analysis.R
#
# Latent Variable Coverage Analysis for the 9-Node because model
# Designed for non-interactive server execution (Rscript coverage_analysis.R)
#
# Scenarios: Varying N_Species, N_Site, N_Surveys
################################################################################

library(because)
library(because.phybase)
library(ape)
library(MASS)
library(parallel)
library(future.apply)

# ---- Configuration -----------------------------------------------------------
N_REPLICATIONS <- 50  # Number of replications per scenario
N_ITER_MCMC <- 12500  # Updated to user-requested iterations
N_BURNIN_MCMC <- 2500 # Adjusted burn-in proportionally
N_CORES_TOTAL <- detectCores() - 1 # Auto-detect server cores
LOG_FILE <- "coverage_9node.log"
RESULTS_FILE <- "coverage_9node_results.rds"

# Ground Truth Parameters for Recovery
TRUE_PARAMS <- list(
  beta_NDVI_U_Resource = 0.5,
  beta_Flower_Cover_U_Resource = 0.6,
  beta_Abundance_Flower_Cover = 0.4,
  beta_Abundance_Temperature_minus_Thermal_Tol = -0.5
)

# ---- Helper: Logging ---------------------------------------------------------
log_msg <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = " "), "\n")
    cat(msg)
    cat(msg, file = LOG_FILE, append = TRUE)
}

# ---- Data Generating Function ------------------------------------------------
generate_9node_data <- function(n_species = 50, n_sites = 30, n_surveys_per_site = 3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # 1. Phylogenetic Tree
  tree <- rtree(n_species)
  tree$tip.label <- paste0("Sp_", 1:n_species)
  
  # 2. Site/Survey structure
  n_surveys <- n_sites * n_surveys_per_site
  sites <- paste0("Site_", 1:n_sites)
  surveys <- paste0("Surv_", 1:n_surveys)
  
  # Environmental variables (Site level)
  Elevation_s <- rnorm(n_sites)
  Temperature <- -1.5 * Elevation_s + rnorm(n_sites, 0, 0.5)
  
  # Latent Resource (Site level)
  U_Resource <- rnorm(n_sites)
  
  # Site-level indicators of U_Resource
  NDVI <- 0.5 * U_Resource + rnorm(n_sites, 0, 0.4)
  Flower_Cover <- 0.6 * U_Resource + 0.4 * Elevation_s + rnorm(n_sites, 0, 0.4)
  
  # 3. Species Traits
  Body_Mass_s <- rnorm(n_species)
  Metabolic_Rate <- 0.8 * Body_Mass_s + rnorm(n_species, 0, 0.3)
  Thermal_Tol <- 0.6 * Body_Mass_s + rnorm(n_species, 0, 0.3)
  
  # 4. Survey level
  Wind_Speed <- rnorm(n_surveys)
  site_idx <- rep(1:n_sites, each = n_surveys_per_site)
  
  # 5. Observations (Full Cross)
  obs_data <- expand.grid(Survey = surveys, Species = tree$tip.label, stringsAsFactors = FALSE)
  
  # Data frames for merging
  d_site <- data.frame(Site = sites, Elevation_s, Temperature, U_Resource, NDVI, Flower_Cover, stringsAsFactors = FALSE)
  d_survey <- data.frame(Survey = surveys, Site = sites[site_idx], Wind_Speed, stringsAsFactors = FALSE)
  d_species <- data.frame(Species = tree$tip.label, Body_Mass_s, Metabolic_Rate, Thermal_Tol, stringsAsFactors = FALSE)
  
  # Merging in sequence
  obs_data <- merge(obs_data, d_survey, by="Survey")
  obs_data <- merge(obs_data, d_species, by="Species")
  obs_data <- merge(obs_data, d_site, by="Site")
  
  # Derived term
  obs_data$Temp_minus_Thermal <- obs_data$Temperature - obs_data$Thermal_Tol
  
  # Response: Abundance (Poisson)
  # beta_Abundance_Flower_Cover = 0.4
  # beta_Abundance_Temp_minus_Thermal = -0.5
  # beta_Abundance_Wind = -0.2
  log_lambda <- -0.5 + 
                0.4 * obs_data$Flower_Cover - 
                0.5 * obs_data$Temp_minus_Thermal - 
                0.2 * obs_data$Wind_Speed +
                rnorm(nrow(obs_data), 0, 0.2) # residual noise
  
  obs_data$Abundance <- rpois(nrow(obs_data), exp(log_lambda))
  
  return(list(
    data = obs_data,
    tree = tree,
    species = d_species,
    site = d_site,
    survey = d_survey
  ))
}

# ---- Simulation Logic --------------------------------------------------------
run_replication <- function(rep_id, n_sp, n_si, n_su) {
  log_msg(sprintf("Starting Rep %d (Sp=%d, Si=%d, Su=%d)", rep_id, n_sp, n_si, n_su))
  
  sim <- generate_9node_data(n_species = n_sp, n_sites = n_si, n_surveys_per_site = n_su, seed = 1000 + rep_id)
  
  eqs <- list(
    NDVI ~ U_Resource,
    Flower_Cover ~ U_Resource + Elevation_s,
    Metabolic_Rate ~ Body_Mass_s,
    Thermal_Tol ~ Body_Mass_s,
    Temperature ~ Elevation_s,
    Abundance ~ Flower_Cover + Wind_Speed + I(Temperature - Thermal_Tol)
  )
  
  t_start <- proc.time()["elapsed"]
  fit <- tryCatch({
    because(
      eqs,
      data = sim$data,
      latent = "U_Resource",
      family = c(Abundance = "poisson"),
      structure = sim$tree,
      n.iter = N_ITER_MCMC,
      n.burnin = N_BURNIN_MCMC,
      quiet = TRUE,
      n.chains = 3,
      n.cores = 1 # Run chains sequentially within worker to avoid nested parallelism issues
    )
  }, error = function(e) {
    log_msg("REPLICATION FAILED:", e$message)
    return(NULL)
  })
  t_end <- proc.time()["elapsed"]
  
  if (is.null(fit)) return(NULL)
  
  summ <- summary(fit)
  res <- summ$results
  
  # Extract focal parameters
  out <- data.frame(
    rep_id = rep_id,
    n_sp = n_sp, n_si = n_si, n_su = n_su,
    time = t_end - t_start
  )
  
  for (p in names(TRUE_PARAMS)) {
    if (p %in% rownames(res)) {
      post_mean <- res[p, "Mean"]
      post_lo <- res[p, "2.5%"]
      post_hi <- res[p, "97.5%"]
      post_rhat <- res[p, "Rhat"]
      
      out[[paste0(p, "_mean")]] <- post_mean
      out[[paste0(p, "_lo")]] <- post_lo
      out[[paste0(p, "_hi")]] <- post_hi
      out[[paste0(p, "_rhat")]] <- post_rhat
      out[[paste0(p, "_covered")]] <- as.integer(TRUE_PARAMS[[p]] >= post_lo && TRUE_PARAMS[[p]] <= post_hi)
      out[[paste0(p, "_bias")]] <- post_mean - TRUE_PARAMS[[p]]
    }
  }
  
  return(out)
}

# ---- Scenario Grid -----------------------------------------------------------
# Middle ground: n_species=50, n_sites=30, n_surveys=3
scenarios <- expand.grid(
  n_species = c(25, 50, 100, 200),
  n_sites = c(15, 30, 60),
  n_surveys = 3,
  replicate = 1:N_REPLICATIONS
)

log_msg("Starting Simulation with scenarios (total replications =", nrow(scenarios), "):")
print(head(scenarios))

# ---- Execution ---------------------------------------------------------------
plan(multisession, workers = N_CORES_TOTAL)

results_list <- future_lapply(1:nrow(scenarios), function(i) {
  scen <- scenarios[i, ]
  run_replication(scen$replicate, scen$n_species, scen$n_sites, scen$n_surveys)
}, future.seed = TRUE)

# Filter out failures
results_df <- do.call(rbind, results_list)

# ---- Save Results -------------------------------------------------------------
saveRDS(results_df, RESULTS_FILE)
log_msg("Simulation Complete. Results saved to", RESULTS_FILE)

# Simple Summary
if (!is.null(results_df)) {
  summary_stats <- aggregate(cbind(beta_NDVI_U_Resource_covered, beta_NDVI_U_Resource_bias) ~ n_species + n_sites, 
                             data = results_df, 
                             FUN = function(x) c(mean = mean(x, na.rm=TRUE)))
  print(summary_stats)
}
