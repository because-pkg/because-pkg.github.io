library(because)
library(because.phybase)
devtools::load_all(
    "/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because"
)
devtools::load_all(
    "/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.phybase"
)

data(rhino.dat)
data(rhino.tree)

sem8_eq <- list(
    LS ~ BM,
    NL ~ BM + RS,
    DD ~ NL
)

message("\n======================================")
message("Fitting model 8 with JAGS (12,000 iterations, 3 parallel chains)...")
message("======================================")
time_jags <- system.time({
    fit_jags <- because(
        equations = sem8_eq,
        data = rhino.dat,
        structure = rhino.tree,
        id_col = "SP",
        engine = 'jags',
        parallel = TRUE,
        n.cores = 3,
        n.chains = 3,
        n.iter = 12000,
        n.burnin = 2000,
        WAIC = FALSE,
        quiet = TRUE
    )
})

message("\n======================================")
message(
    "Fitting model 8 with NIMBLE + AF_slice (12,000 iterations, 3 parallel chains)..."
)
message("AF_slice matches JAGS accuracy; parallelization matches JAGS speed.")
message("======================================")
time_nimble <- system.time({
    fit_nimble <- because(
        equations = sem8_eq,
        data = rhino.dat,
        structure = rhino.tree,
        id_col = "SP",
        engine = 'nimble',
        parallel = TRUE,
        n.cores = 3,
        n.chains = 3,
        n.iter = 12000,
        n.burnin = 2000,
        WAIC = FALSE,
        quiet = FALSE
    )
})

message("\n======================================")
message("RESULTS COMPARISON: Rhinograd Model 8")
message("(JAGS: 12k iters PARALLEL | NIMBLE+AF_slice: 12k iters PARALLEL)")
message("======================================")

params <- rownames(fit_jags$summary$statistics)

cat("\n--- PARAMETER ESTIMATES (Mean) ---\n")
means_df <- data.frame(
    JAGS_12k = fit_jags$summary$statistics[params, "Mean"],
    NIMBLE_3k = fit_nimble$summary$statistics[params, "Mean"]
)
print(means_df)

cat("\n--- PARAMETER ESTIMATES (SD) ---\n")
sd_df <- data.frame(
    JAGS_12k = fit_jags$summary$statistics[params, "SD"],
    NIMBLE_3k = fit_nimble$summary$statistics[params, "SD"]
)
print(sd_df)

cat("\n--- EXECUTION TIME ---\n")
time_df <- data.frame(
    JAGS_12k_sec = time_jags["elapsed"],
    NIMBLE_3k_sec = time_nimble["elapsed"]
)
rownames(time_df) <- "Total Elapsed"
print(time_df)

message("\nDone!")
