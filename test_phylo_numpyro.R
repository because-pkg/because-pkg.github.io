library(because)
library(because.phybase)

data(rhino.dat)
data(rhino.tree)

sem8_eq <- list(
  LS ~ BM,
  NL ~ BM + RS,
  DD ~ NL
)

cat("Fitting JAGS...\n")
fit_sem8 <- because(
  equations = sem8_eq,
  data = rhino.dat,
  structure = rhino.tree,
  id_col = "SP",
  WAIC = TRUE,
  parallel = TRUE,
  n.cores = 3,
  n.iter = 1000,
  n.burnin = 500,
  engine = "jags"
)

cat("Fitting NumPyro...\n")
fit_sem8_py <- because(
  equations = sem8_eq,
  data = rhino.dat,
  structure = rhino.tree,
  id_col = "SP",
  WAIC = TRUE,
  parallel = TRUE,
  n.cores = 3,
  n.iter = 1000,
  n.burnin = 500,
  engine = "numpyro"
)

cat("\n--- JAGS Summary ---\n")
print(summary(fit_sem8))

cat("\n--- NumPyro Summary ---\n")
print(summary(fit_sem8_py))
