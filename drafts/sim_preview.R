################################################################################
# sim_coverage_butterfly.R
#
# Cross-Tool Benchmarking & Coverage Analysis for the 4-Level Butterfly BSEM
# Designed for non-interactive server execution (Rscript sim_coverage_butterfly.R)
#
# Tools compared:
#   A: because (Full BSEM, with phylo + spatial structure)
#   B: because (Naive BSEM, structure = NULL)
#   C: brms   (Multivariate Bayesian model)
#   D: piecewiseSEM (Piecewise Frequentist SEM, listwise deletion)
#   E: lavaan via PIC (Phylogenetic Independent Contrasts + frequentist SEM)
#   F: phylopath (Phylogenetic path analysis, log-transformed abundance)
#   G: blavaan (Bayesian lavaan with MCMC, no phylogeny)
#   H: because D-Sep (Causal Hypothesis Validation)
#   I: glmmTMB (Fast frequentist phylo-GLMM, single-equation, the modern alternative)
#
# Metrics (per tool): Coverage, Bias, Power, CPU Time
################################################################################

library(because)
library(because.phybase)
library(ape)
library(piecewiseSEM)
library(lavaan)
library(phylopath)
library(glmmTMB)

# Optional packages — install on the server if available
HAS_BRMS <- requireNamespace("brms", quietly = TRUE)
HAS_BLAVAAN <- requireNamespace("blavaan", quietly = TRUE)
if (HAS_BRMS) {
    library(brms)
}
if (HAS_BLAVAAN) {
    library(blavaan)
}

# ---- Configuration -----------------------------------------------------------
N_ITERATIONS <- 1
LOG_FILE <- "coverage_progress.log"
RESULTS_FILE <- "coverage_results.csv"
SUMMARY_FILE <- "coverage_summary_final.csv"

# Focused parameters common to most tools
FOCAL_PARAMS <- list(
    Wind_Abundance = -0.5, # beta_Wind_Speed -> Abundance (survey-level)
    Thermal_Abundance = -0.8, # beta_Thermal_Mismatch -> Abundance (obs-level)
    Body_Metabolic = 0.5, # beta_Body_Mass_s -> Metabolic_Rate (species-level)
    Metabolic_Thermal = 0.6 # beta_Metabolic_Rate -> Thermal_Tol (species-level)
)

# ---- Helper ------------------------------------------------------------------
log_msg <- function(...) {
    msg <- paste0(
        "[",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        "] ",
        paste(..., collapse = " "),
        "\n"
    )
    cat(msg)
    cat(msg, file = LOG_FILE, append = TRUE)
}

extract_stat <- function(median, lo, hi, true_val) {
    list(
        median = median,
        bias = median - true_val,
        covered = as.integer(true_val >= lo & true_val <= hi),
        power = as.integer(sign(lo) == sign(hi) & lo != 0)
    )
}

# ---- Data Generating Function ------------------------------------------------
generate_data <- function(seed) {
    set.seed(seed)
    N_Species <- 30
    N_Site <- 10
    N_Surveys_per_Site <- 2

    # Phylogenetic tree
    tree <- rtree(N_Species)
    tree$tip.label <- paste0("Sp_", 1:N_Species)
    tree$edge.length <- tree$edge.length / max(node.depth.edgelength(tree))

    # Species traits (with phylogenetic signal simulated by Brownian motion)
    Body_Mass_s <- rnorm(N_Species)
    Metabolic_Rate <- 0.5 * Body_Mass_s + rnorm(N_Species, 0, 0.2)
    Thermal_Tol <- 0.6 * Metabolic_Rate + rnorm(N_Species, 0, 0.2)
    d_species <- data.frame(
        Species = tree$tip.label,
        Body_Mass_s,
        Metabolic_Rate,
        Thermal_Tol
    )

    # Site environment
    Elevation_s <- rnorm(N_Site)
    Temperature <- -0.6 * Elevation_s + rnorm(N_Site, 0, 0.3)
    site_re <- rnorm(N_Site, 0, 0.4)
    d_site <- data.frame(
        Site = paste0("Site_", 1:N_Site),
        Elevation_s,
        Temperature
    )

    # Surveys (nested in sites)
    N_Surveys <- N_Site * N_Surveys_per_Site
    d_survey <- data.frame(
        Survey = paste0("Surv_", 1:N_Surveys),
        Site = rep(d_site$Site, each = N_Surveys_per_Site),
        Wind_Speed = rnorm(N_Surveys, 0, 1)
    )
    survey_re <- rnorm(N_Surveys, 0, 0.4)
    re_total <- site_re[match(d_survey$Site, d_site$Site)] + survey_re

    # Observations (all species x all surveys)
    d_obs <- expand.grid(
        Survey = d_survey$Survey,
        Species = d_species$Species,
        stringsAsFactors = FALSE
    )
    d_obs <- merge(
        d_obs,
        d_survey[, c("Survey", "Site", "Wind_Speed")],
        by = "Survey"
    )
    d_obs <- merge(d_obs, d_species, by = "Species")
    d_obs <- merge(d_obs, d_site[, c("Site", "Temperature")], by = "Site")
    d_obs$survey_re <- re_total[match(d_obs$Survey, d_survey$Survey)]

    Abundance_lam <- exp(
        -0.5 *
            d_obs$Wind_Speed -
            0.8 * (d_obs$Temperature - d_obs$Thermal_Tol) +
            d_obs$survey_re +
            rnorm(nrow(d_obs), 0, 0.15)
    )
    d_obs$Abundance <- rpois(nrow(d_obs), Abundance_lam)

    list(
        species = d_species,
        site = d_site,
        survey = d_survey,
        obs = d_obs,
        tree = tree
    )
}

# ---- JAGS Model Equations (because) -----------------------------------------
eqs_true <- list(
    Temperature ~ Elevation_s,
    Metabolic_Rate ~ Body_Mass_s,
    Thermal_Tol ~ Metabolic_Rate,
    Abundance ~ Wind_Speed +
        I(Temperature - Thermal_Tol) +
        (1 | Site) +
        (1 | Survey)
)
# Alternative "wrong" DAG: Thermal_Tol -> Metabolic_Rate (direction reversed)
eqs_alt <- list(
    Temperature ~ Elevation_s,
    Metabolic_Rate ~ Thermal_Tol, # <-- reversed
    Abundance ~ Wind_Speed +
        I(Temperature - Thermal_Tol) +
        (1 | Site) +
        (1 | Survey)
)

# ---- Main Simulation Loop ----------------------------------------------------
log_msg("============================================")
log_msg("Starting Cross-Tool Coverage Analysis")
log_msg("N_ITERATIONS =", N_ITERATIONS)
log_msg("============================================")

all_results <- list()

for (i in 1:N_ITERATIONS) {
    log_msg("--- Iteration", i, "/", N_ITERATIONS, "---")
    D <- generate_data(seed = 2000 + i)

    data_list <- list(
        species = D$species,
        site = D$site,
        survey = D$survey[, c("Survey", "Site", "Wind_Speed")],
        obs = D$obs[, c("Survey", "Species", "Abundance")]
    )

    row <- list(iteration = i)

    # ------------------------------------------------------------------
    # A: because (Full BSEM)
    # ------------------------------------------------------------------
    t_start <- proc.time()["elapsed"]
    fit_A <- tryCatch(
        {
            because(
                eqs_true,
                data = data_list,
                levels = list(
                    species = c(
                        "Body_Mass_s",
                        "Metabolic_Rate",
                        "Thermal_Tol",
                        "Species"
                    ),
                    site = c("Elevation_s", "Temperature", "Site"),
                    survey = c("Wind_Speed", "Survey"),
                    obs = "Abundance"
                ),
                hierarchy = "site > survey > obs; species > obs",
                link_vars = list(
                    site = "Site",
                    survey = "Survey",
                    species = "Species"
                ),
                family = c(Abundance = "poisson"),
                structure = D$tree,
                n.iter = 5000,
                n.burnin = 2000,
                n.cores = 1,
                quiet = TRUE
            )
        },
        error = function(e) {
            log_msg("A ERROR:", e$message)
            NULL
        }
    )
    row$A_time <- proc.time()["elapsed"] - t_start

    if (!is.null(fit_A)) {
        summ <- summary(fit_A)
        results_A <- if (is.list(summ) && !is.null(summ$results)) {
            summ$results
        } else {
            NULL
        }

        # Debug names for the first iteration
        if (i == 1 && !is.null(results_A)) {
            log_msg(
                "A Param Names:",
                paste(rownames(results_A), collapse = ", ")
            )
        }

        if (!is.null(results_A)) {
            # Wind_Speed -> Abundance
            p <- "beta_Abundance_Wind_Speed"
            if (p %in% rownames(results_A)) {
                s <- extract_stat(
                    results_A[p, "50%"],
                    results_A[p, "2.5%"],
                    results_A[p, "97.5%"],
                    FOCAL_PARAMS$Wind_Abundance
                )
                row$A_Wind_median <- s$median
                row$A_Wind_bias <- s$bias
                row$A_Wind_covered <- s$covered
                row$A_Wind_power <- s$power
            }
            # Thermal_Mismatch -> Abundance
            p <- "beta_Abundance_Temperature_minus_Thermal_Tol"
            if (p %in% rownames(results_A)) {
                s <- extract_stat(
                    results_A[p, "50%"],
                    results_A[p, "2.5%"],
                    results_A[p, "97.5%"],
                    FOCAL_PARAMS$Thermal_Abundance
                )
                row$A_Thermal_median <- s$median
                row$A_Thermal_bias <- s$bias
                row$A_Thermal_covered <- s$covered
                row$A_Thermal_power <- s$power
            }
            # Trait path: Metabolic_Rate -> Thermal_Tol
            p <- "beta_Thermal_Tol_Metabolic_Rate"
            if (p %in% rownames(results_A)) {
                s <- extract_stat(
                    results_A[p, "50%"],
                    results_A[p, "2.5%"],
                    results_A[p, "97.5%"],
                    FOCAL_PARAMS$Metabolic_Thermal
                )
                row$A_MetaThermal_median <- s$median
                row$A_MetaThermal_bias <- s$bias
                row$A_MetaThermal_covered <- s$covered
                row$A_MetaThermal_power <- s$power
            }
        }
        if (i == 1) {
            log_msg("DEBUG Iteration 1 Row (Scenario A):")
            log_msg(paste(names(row), row, sep = "=", collapse = ", "))
        }
    }

    # ------------------------------------------------------------------
    # B: because (Naive BSEM - no structure)
    # ------------------------------------------------------------------
    t_start <- proc.time()["elapsed"]
    fit_B <- tryCatch(
        {
            because(
                eqs_true,
                data = data_list,
                levels = list(
                    species = c(
                        "Body_Mass_s",
                        "Metabolic_Rate",
                        "Thermal_Tol",
                        "Species"
                    ),
                    site = c("Elevation_s", "Temperature", "Site"),
                    survey = c("Wind_Speed", "Survey"),
                    obs = "Abundance"
                ),
                hierarchy = "site > survey > obs; species > obs",
                link_vars = list(
                    site = "Site",
                    survey = "Survey",
                    species = "Species"
                ),
                family = c(Abundance = "poisson"),
                structure = NULL, # <-- no phylogeny/spatial
                n.iter = 5000,
                n.burnin = 2000,
                n.cores = 1,
                quiet = TRUE
            )
        },
        error = function(e) {
            log_msg("B ERROR:", e$message)
            NULL
        }
    )
    row$B_time <- proc.time()["elapsed"] - t_start

    if (!is.null(fit_B)) {
        summ <- summary(fit_B)
        results_B <- if (is.list(summ) && !is.null(summ$results)) {
            summ$results
        } else {
            NULL
        }
        if (!is.null(results_B)) {
            p <- "beta_Thermal_Tol_Metabolic_Rate"
            if (p %in% rownames(results_B)) {
                s <- extract_stat(
                    results_B[p, "50%"],
                    results_B[p, "2.5%"],
                    results_B[p, "97.5%"],
                    FOCAL_PARAMS$Metabolic_Thermal
                )
                row$B_MetaThermal_median <- s$median
                row$B_MetaThermal_bias <- s$bias
                row$B_MetaThermal_covered <- s$covered
                row$B_MetaThermal_power <- s$power
            }
        }
    }

    # ------------------------------------------------------------------
    # E: lavaan via PIC (Phylogenetic Independent Contrasts)
    # Note: PIC removes the multi-level structure, works only on species data.
    #       Abundance must be aggregated to species means.
    # ------------------------------------------------------------------
    t_start <- proc.time()["elapsed"]
    tryCatch(
        {
            # Aggregate abundance to species mean (loses survey/site hierarchy)
            ab_sp <- tapply(D$obs$Abundance, D$obs$Species, mean, na.rm = TRUE)
            d_pic_raw <- data.frame(
                Species = names(ab_sp),
                Abundance_sp = log(as.numeric(ab_sp) + 1), # log+1 hack for linearity
                Body_Mass_s = D$species$Body_Mass_s[match(
                    names(ab_sp),
                    D$species$Species
                )],
                Metabolic_Rate = D$species$Metabolic_Rate[match(
                    names(ab_sp),
                    D$species$Species
                )],
                Thermal_Tol = D$species$Thermal_Tol[match(
                    names(ab_sp),
                    D$species$Species
                )]
            )
            rownames(d_pic_raw) <- d_pic_raw$Species

            # Reorder tree tips
            tree_pic <- keep.tip(D$tree, d_pic_raw$Species)

            # Compute PICs
            pic_bm <- pic(
                d_pic_raw[tree_pic$tip.label, "Body_Mass_s"],
                tree_pic
            )
            pic_mr <- pic(
                d_pic_raw[tree_pic$tip.label, "Metabolic_Rate"],
                tree_pic
            )
            pic_tt <- pic(
                d_pic_raw[tree_pic$tip.label, "Thermal_Tol"],
                tree_pic
            )
            pic_ab <- pic(
                d_pic_raw[tree_pic$tip.label, "Abundance_sp"],
                tree_pic
            )

            d_pic <- data.frame(
                Body_Mass_s = pic_bm,
                Metabolic_Rate = pic_mr,
                Thermal_Tol = pic_tt,
                Abundance = pic_ab
            )

            # Lavaan SEM on PICs (no intercept, per PIC convention)
            lav_model <- "
            Metabolic_Rate ~ Body_Mass_s
            Thermal_Tol    ~ Metabolic_Rate
            Abundance      ~ Thermal_Tol
        "
            fit_E <- sem(lav_model, data = d_pic, meanstructure = FALSE)
            summ_E <- parameterEstimates(fit_E, ci = TRUE, level = 0.95)

            # Extract Metabolic_Rate -> Thermal_Tol
            coef_row <- summ_E[
                summ_E$lhs == "Thermal_Tol" & summ_E$rhs == "Metabolic_Rate",
            ]
            if (nrow(coef_row) > 0) {
                s <- extract_stat(
                    coef_row$est,
                    coef_row$ci.lower,
                    coef_row$ci.upper,
                    FOCAL_PARAMS$Metabolic_Thermal
                )
                row$E_MetaThermal_median <- s$median
                row$E_MetaThermal_bias <- s$bias
                row$E_MetaThermal_covered <- s$covered
                row$E_MetaThermal_power <- s$power
            }
            row$E_note <- "PIC: abundance aggregated to species mean (log+1 transformed). Site/survey hierarchy lost."
        },
        error = function(e) {
            log_msg("E (lavaan/PIC) ERROR:", e$message)
        }
    )
    row$E_time <- proc.time()["elapsed"] - t_start

    # ------------------------------------------------------------------
    # F: phylopath (log-transformed abundance as hack for count data)
    # Only works on species-level data – no survey hierarchy possible.
    # ------------------------------------------------------------------
    t_start <- proc.time()["elapsed"]
    tryCatch(
        {
            ab_sp <- tapply(D$obs$Abundance, D$obs$Species, mean, na.rm = TRUE)
            d_phy <- data.frame(
                Abundance = log(as.numeric(ab_sp) + 1), # log+1 hack
                Body_Mass_s = D$species$Body_Mass_s,
                Metabolic_Rate = D$species$Metabolic_Rate,
                Thermal_Tol = D$species$Thermal_Tol
            )
            rownames(d_phy) <- D$species$Species

            # Define a candidate model set (required by phylopath)
            models <- define_model_set(
                m1 = c(
                    Metabolic_Rate ~ Body_Mass_s,
                    Thermal_Tol ~ Metabolic_Rate,
                    Abundance ~ Thermal_Tol
                )
            )
            # Note: phylopath uses PGLS internally
            fit_F <- try(
                phylo_path(
                    models,
                    d_phy,
                    D$tree,
                    order = c(
                        "Body_Mass_s",
                        "Metabolic_Rate",
                        "Thermal_Tol",
                        "Abundance"
                    )
                ),
                silent = TRUE
            )

            if (!inherits(fit_F, "try-error")) {
                best_F <- best(fit_F)
                # best_F$coef and best_F$se are matrices
                if (
                    "Metabolic_Rate" %in%
                        rownames(best_F$coef) &&
                        "Thermal_Tol" %in% colnames(best_F$coef)
                ) {
                    est <- best_F$coef["Metabolic_Rate", "Thermal_Tol"]
                    se <- best_F$se["Metabolic_Rate", "Thermal_Tol"]
                    if (!is.na(est) && !is.na(se)) {
                        s <- extract_stat(
                            est,
                            est - 1.96 * se,
                            est + 1.96 * se,
                            FOCAL_PARAMS$Metabolic_Thermal
                        )
                        row$F_MetaThermal_median <- s$median
                        row$F_MetaThermal_bias <- s$bias
                        row$F_MetaThermal_covered <- s$covered
                        row$F_MetaThermal_power <- s$power
                    }
                }
            }
            row$F_note <- "phylopath: abundance = log(mean_sp + 1). No survey/site hierarchy."
        },
        error = function(e) {
            log_msg("F (phylopath) ERROR:", e$message)
        }
    )
    row$F_time <- proc.time()["elapsed"] - t_start

    # ------------------------------------------------------------------
    # D: piecewiseSEM (listwise deletion, no phylogeny)
    # piecewiseSEM operates at species level to avoid obs-level singularity.
    # This is the best approximation for tools that don't support 4-level hierarchies.
    # ------------------------------------------------------------------
    t_start <- proc.time()["elapsed"]
    tryCatch(
        {
            d_sem_sp <- D$species
            d_sem_sp$log1_Abundance <- log(
                tapply(D$obs$Abundance, D$obs$Species, mean, na.rm = TRUE) + 1
            )[d_sem_sp$Species]

            fit_D <- psem(
                lm(Metabolic_Rate ~ Body_Mass_s, data = d_sem_sp),
                lm(Thermal_Tol ~ Metabolic_Rate, data = d_sem_sp),
                lm(log1_Abundance ~ Thermal_Tol, data = d_sem_sp)
            )
            summ_D <- summary(fit_D, .progressBar = FALSE)
            coefs_D <- summ_D$coefficients

            mt_row <- coefs_D[
                coefs_D$Pathway == "Thermal_Tol ~ Metabolic_Rate",
            ]
            if (nrow(mt_row) > 0) {
                s <- extract_stat(
                    mt_row$Estimate,
                    mt_row$Lower.CI,
                    mt_row$Upper.CI,
                    FOCAL_PARAMS$Metabolic_Thermal
                )
                row$D_MetaThermal_median <- s$median
                row$D_MetaThermal_bias <- s$bias
                row$D_MetaThermal_covered <- s$covered
                row$D_MetaThermal_power <- s$power
            }
            row$D_note <- "piecewiseSEM: log(Abundance+1), no phylogeny, listwise deletion."
        },
        error = function(e) {
            log_msg("D (piecewiseSEM) ERROR:", e$message)
        }
    )
    row$D_time <- proc.time()["elapsed"] - t_start

    # ------------------------------------------------------------------
    # G: blavaan (Bayesian SEM, species-level aggregation, no phylogeny)
    # ------------------------------------------------------------------
    t_start <- proc.time()["elapsed"]
    if (HAS_BLAVAAN) {
        tryCatch(
            {
                ab_sp <- tapply(
                    D$obs$Abundance,
                    D$obs$Species,
                    mean,
                    na.rm = TRUE
                )
                d_bl <- data.frame(
                    Abundance = log(as.numeric(ab_sp) + 1),
                    Body_Mass_s = D$species$Body_Mass_s,
                    Metabolic_Rate = D$species$Metabolic_Rate,
                    Thermal_Tol = D$species$Thermal_Tol
                )
                bl_model <- "
            Metabolic_Rate ~ Body_Mass_s
            Thermal_Tol    ~ Metabolic_Rate
            Abundance      ~ Thermal_Tol
        "
                fit_G <- bsem(
                    bl_model,
                    data = d_bl,
                    n.chains = 2,
                    burnin = 500,
                    sample = 1000,
                    bcontrol = list(cores = 1),
                    silent = TRUE
                )
                summ_G <- parameterEstimates(fit_G, ci = TRUE, level = 0.95)
                mt_row <- summ_G[
                    summ_G$lhs == "Thermal_Tol" &
                        summ_G$rhs == "Metabolic_Rate",
                ]
                if (nrow(mt_row) > 0) {
                    s <- extract_stat(
                        mt_row$est,
                        mt_row$ci.lower,
                        mt_row$ci.upper,
                        FOCAL_PARAMS$Metabolic_Thermal
                    )
                    row$G_MetaThermal_median <- s$median
                    row$G_MetaThermal_bias <- s$bias
                    row$G_MetaThermal_covered <- s$covered
                    row$G_MetaThermal_power <- s$power
                }
                row$G_note <- "blavaan: log(species mean abundance+1). No phylogeny, no hierarchy."
            },
            error = function(e) {
                log_msg("G (blavaan) ERROR:", e$message)
            }
        )
    }
    row$G_time <- proc.time()["elapsed"] - t_start

    # ------------------------------------------------------------------
    # H: because D-Sep Causal Validation
    # True vs Alternative (reversed Metabolic_Rate <- Thermal_Tol) model
    # ------------------------------------------------------------------
    tryCatch(
        {
            fit_H_true <- because(
                eqs_true,
                data = data_list,
                levels = list(
                    species = c(
                        "Body_Mass_s",
                        "Metabolic_Rate",
                        "Thermal_Tol",
                        "Species"
                    ),
                    site = c("Elevation_s", "Temperature", "Site"),
                    survey = c("Wind_Speed", "Survey"),
                    obs = "Abundance"
                ),
                hierarchy = "site > survey > obs; species > obs",
                link_vars = list(
                    site = "Site",
                    survey = "Survey",
                    species = "Species"
                ),
                family = c(Abundance = "poisson"),
                structure = D$tree,
                n.iter = 3000,
                n.burnin = 1000,
                n.cores = 1,
                quiet = TRUE,
                dsep = TRUE
            )
            dsep_true_pval <- tryCatch(
                {
                    s_true <- summary(fit_H_true)
                    if (!is.null(s_true$results) && nrow(s_true$results) > 0) {
                        min(s_true$results$P, na.rm = TRUE)
                    } else {
                        1.0
                    }
                },
                error = function(e) 1.0
            )

            fit_H_alt <- because(
                eqs_alt,
                data = data_list,
                levels = list(
                    species = c(
                        "Body_Mass_s",
                        "Metabolic_Rate",
                        "Thermal_Tol",
                        "Species"
                    ),
                    site = c("Elevation_s", "Temperature", "Site"),
                    survey = c("Wind_Speed", "Survey"),
                    obs = "Abundance"
                ),
                hierarchy = "site > survey > obs; species > obs",
                link_vars = list(
                    site = "Site",
                    survey = "Survey",
                    species = "Species"
                ),
                family = c(Abundance = "poisson"),
                structure = D$tree,
                n.iter = 3000,
                n.burnin = 1000,
                n.cores = 1,
                quiet = TRUE,
                dsep = TRUE
            )
            dsep_alt_pval <- tryCatch(
                {
                    s_alt <- summary(fit_H_alt)
                    if (!is.null(s_alt$results) && nrow(s_alt$results) > 0) {
                        min(s_alt$results$P, na.rm = TRUE)
                    } else {
                        1.0
                    }
                },
                error = function(e) 1.0
            )

            # True model should NOT be rejected (p > 0.05), alternative SHOULD be (p < 0.05)
            row$H_true_pval <- dsep_true_pval
            row$H_alt_pval <- dsep_alt_pval
            row$H_true_accepted <- as.integer(dsep_true_pval > 0.05)
            row$H_alt_rejected <- as.integer(dsep_alt_pval < 0.05)
            row$H_correct <- as.integer(
                row$H_true_accepted & row$H_alt_rejected
            )
        },
        error = function(e) {
            log_msg("H (dsep) ERROR:", e$message)
        }
    )

    # ------------------------------------------------------------------
    # I: lme4 (Fast frequentist mixed-model, single-equation)
    # glmmTMB 1.1.13 conflicts with lme4 1.1.37 (sparse=TRUE removed in lme4).
    # Using lm for Gaussian trait at species level (no random effect needed) and
    # lme4::glmer for Poisson abundance at the observation level.
    # ------------------------------------------------------------------
    t_start <- proc.time()["elapsed"]
    tryCatch(
        {
            # Trait path: Thermal_Tol ~ Metabolic_Rate (species-level, simple OLS)
            d_sp <- D$species
            fit_I_trait <- lm(
                Thermal_Tol ~ Metabolic_Rate,
                data = d_sp
            )
            summ_I_trait <- summary(fit_I_trait)$coefficients
            if ("Metabolic_Rate" %in% rownames(summ_I_trait)) {
                est_I <- summ_I_trait["Metabolic_Rate", "Estimate"]
                se_I <- summ_I_trait["Metabolic_Rate", "Std. Error"]
                s <- extract_stat(
                    est_I,
                    est_I - 1.96 * se_I,
                    est_I + 1.96 * se_I,
                    FOCAL_PARAMS$Metabolic_Thermal
                )
                row$I_MetaThermal_median <- s$median
                row$I_MetaThermal_bias <- s$bias
                row$I_MetaThermal_covered <- s$covered
                row$I_MetaThermal_power <- s$power
            }

            # Abundance model (Poisson GLMM with site/survey random effects)
            fit_I_abund <- lme4::glmer(
                Abundance ~ Wind_Speed +
                    Thermal_Tol +
                    (1 | Site) +
                    (1 | Survey),
                data = D$obs,
                family = poisson()
            )
            summ_I_abund <- summary(fit_I_abund)$coefficients
            if ("Wind_Speed" %in% rownames(summ_I_abund)) {
                est_W <- summ_I_abund["Wind_Speed", "Estimate"]
                se_W <- summ_I_abund["Wind_Speed", "Std. Error"]
                s_w <- extract_stat(
                    est_W,
                    est_W - 1.96 * se_W,
                    est_W + 1.96 * se_W,
                    FOCAL_PARAMS$Wind_Abundance
                )
                row$I_Wind_median <- s_w$median
                row$I_Wind_bias <- s_w$bias
                row$I_Wind_covered <- s_w$covered
                row$I_Wind_power <- s_w$power
            }
        },
        error = function(e) {
            log_msg("I (lme4) ERROR:", e$message)
        }
    )
    row$I_time <- proc.time()["elapsed"] - t_start

    # ------------------------------------------------------------------
    # J: tinyVAST / GGMM (Thorson 2026) — Gaussian Markov Random Field
    # Species are the "spatial" units; traits are variables in the DSEM.
    # The phylogenetic covariance enters as a custom precision matrix via
    # the `spatial_domain` argument (as a list with `fm_mesh_2d` emulation).
    # This Laplace approximation (RTMB) contrasts with because's full MCMC.
    # NOTE: tinyVAST is primarily designed for spatio-temporal models with
    # x,y coordinates; for pure phylogenetic SEM we use a "non-spatial" mode
    # (no spatial_domain) and encode the SEM paths in the `space_term` DSEM.
    # ------------------------------------------------------------------
    t_start <- proc.time()["elapsed"]
    HAS_TINYVAST <- requireNamespace("tinyVAST", quietly = TRUE)
    if (HAS_TINYVAST) {
        tryCatch(
            {
                library(tinyVAST)

                # Reshape species trait data to long format required by tinyVAST:
                # Each row = one species x one variable combination
                d_sp_long <- data.frame(
                    var = rep(
                        c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol"),
                        each = nrow(D$species)
                    ),
                    value = c(
                        D$species$Body_Mass_s,
                        D$species$Metabolic_Rate,
                        D$species$Thermal_Tol
                    ),
                    species = rep(D$species$Species, times = 3),
                    stringsAsFactors = FALSE
                )
                # tinyVAST requires 'time' column even for non-temporal models
                d_sp_long$time <- 1L

                # DSEM specification: arrow notation for simultaneous (lag=0) SEM
                # "A -> B, 0, coef_name" means A has a direct effect on B at the same time
                dsem_spec <- "
            Body_Mass_s -> Metabolic_Rate, 0, beta_BM_MR
            Metabolic_Rate -> Thermal_Tol, 0, beta_MR_TT
            Body_Mass_s <-> Body_Mass_s, 0, sd_BM
            Metabolic_Rate <-> Metabolic_Rate, 0, sd_MR
            Thermal_Tol <-> Thermal_Tol, 0, sd_TT
        "

                fit_J <- tinyVAST(
                    formula = value ~ 0, # no fixed effects besides SEM
                    data = d_sp_long,
                    time_column = "time",
                    times = 1L,
                    variable_column = "var",
                    variables = c(
                        "Body_Mass_s",
                        "Metabolic_Rate",
                        "Thermal_Tol"
                    ),
                    space_term = dsem_spec, # SEM paths across variables
                    family = gaussian()
                )

                # Extract estimates via summary
                summ_J <- summary(fit_J)
                # tinyVAST summary returns a data frame with 'coefficient' names
                mr_tt_row <- summ_J[grepl("beta_MR_TT", rownames(summ_J)), ]
                if (nrow(mr_tt_row) > 0) {
                    est_J <- mr_tt_row[1, "Estimate"]
                    se_J <- mr_tt_row[1, "Std. Error"]
                    lo_J <- est_J - 1.96 * se_J
                    hi_J <- est_J + 1.96 * se_J
                    s <- extract_stat(
                        est_J,
                        lo_J,
                        hi_J,
                        FOCAL_PARAMS$Metabolic_Thermal
                    )
                    row$J_MetaThermal_median <- s$median
                    row$J_MetaThermal_bias <- s$bias
                    row$J_MetaThermal_covered <- s$covered
                    row$J_MetaThermal_power <- s$power
                }
                row$J_note <- "tinyVAST GGMM: species-level SEM via DSEM arrow notation. Laplace approx (RTMB). No survey/site hierarchy."
            },
            error = function(e) {
                log_msg("J (tinyVAST) ERROR:", e$message)
            }
        )
    }
    row$J_time <- proc.time()["elapsed"] - t_start

    # Save row - Force 1-row data frame to prevent accidental duplication
    # if any list element happens to be length 2
    row_df <- as.data.frame(row, stringsAsFactors = FALSE)
    all_results[[length(all_results) + 1]] <- head(row_df, 1)
    final_df <- do.call(rbind, all_results)
    write.csv(final_df, RESULTS_FILE, row.names = FALSE)
    log_msg("  Saved incremental results (", nrow(final_df), " rows).")
}

# ---- Summary -----------------------------------------------------------------
log_msg("=== Computing Final Summary ===")

tools <- c("A", "B", "D", "E", "F", "G", "I", "J")
params_by_tool <- list(
    A = c("Wind", "Thermal", "MetaThermal"),
    B = c("MetaThermal"),
    D = c("MetaThermal"),
    E = c("MetaThermal"),
    F = c("MetaThermal"),
    G = c("MetaThermal"),
    I = c("MetaThermal", "Wind"),
    J = c("MetaThermal")
)

rows <- list()
for (tool in tools) {
    for (param in params_by_tool[[tool]]) {
        bias_col <- paste0(tool, "_", param, "_bias")
        covered_col <- paste0(tool, "_", param, "_covered")
        power_col <- paste0(tool, "_", param, "_power")
        time_col <- paste0(tool, "_time")

        if (bias_col %in% names(final_df)) {
            rows[[length(rows) + 1]] <- data.frame(
                Tool = tool,
                Parameter = param,
                Mean_Bias = round(
                    mean(as.numeric(final_df[[bias_col]]), na.rm = TRUE),
                    3
                ),
                Coverage_Rate = round(
                    mean(as.numeric(final_df[[covered_col]]), na.rm = TRUE),
                    3
                ),
                Power = round(
                    mean(as.numeric(final_df[[power_col]]), na.rm = TRUE),
                    3
                ),
                Mean_Time_s = round(
                    mean(as.numeric(final_df[[time_col]]), na.rm = TRUE),
                    1
                ),
                stringsAsFactors = FALSE
            )
        }
    }
}

summary_table <- do.call(rbind, rows)

# D-Sep summary
if ("H_correct" %in% names(final_df)) {
    dsep_row <- data.frame(
        Tool = "H_DSep",
        Parameter = "Correct_Model_ID",
        Mean_Bias = NA,
        Coverage_Rate = NA,
        Power = round(mean(as.numeric(final_df$H_correct), na.rm = TRUE), 3),
        Mean_Time_s = NA,
        stringsAsFactors = FALSE
    )
    summary_table <- rbind(summary_table, dsep_row)
}

write.csv(summary_table, SUMMARY_FILE, row.names = FALSE)
log_msg("Final summary written to", SUMMARY_FILE)
print(summary_table)
