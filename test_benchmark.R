library(because)
library(because.phybase)

data(rhino.dat)
data(rhino.tree)

sem8_eq <- list(
    LS ~ BM,
    NL ~ BM + RS,
    DD ~ NL
)

cat("\n=== RUNNING JAGS ===\n")
t1 <- Sys.time()
fit_jags <- because(
    equations = sem8_eq,
    data = rhino.dat,
    structure = rhino.tree,
    id_col = "SP",
    engine = "jags",
    parallel = TRUE,
    n.cores = 3,
    n.iter = 2000,
    n.burnin = 1000
)
t2 <- Sys.time()
cat(sprintf("JAGS time: %.2f secs\n", as.numeric(difftime(t2, t1, units="secs"))))

cat("\n=== RUNNING NUMPYRO ===\n")
t1 <- Sys.time()
fit_numpyro <- because(
    equations = sem8_eq,
    data = rhino.dat,
    structure = rhino.tree,
    id_col = "SP",
    engine = "numpyro",
    parallel = TRUE,
    n.cores = 3,
    n.iter = 2000,
    n.warmup = 1000
)
t2 <- Sys.time()
cat(sprintf("NumPyro time: %.2f secs\n", as.numeric(difftime(t2, t1, units="secs"))))

cat("\n=== JAGS SUMMARY ===\n")
print(summary(fit_jags))

cat("\n=== NUMPYRO SUMMARY ===\n")
print(summary(fit_numpyro))

