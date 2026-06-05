library(because.phybase)
data(rhino.dat)
data(rhino.tree)
sem8_eq <- list(
  LS ~ BM,
  NL ~ BM + RS,
  DD ~ NL
)

# Use trace to print equations right before python call
trace(because:::because, edit = FALSE, at = list(c(67, 4)), tracer = quote({
  if (engine == "numpyro") {
     cat("=== EQUATIONS PASSED ===\n")
     print(equations)
     cat("========================\n")
  }
}))

fit <- because(equations = sem8_eq, data = rhino.dat, structure = rhino.tree, id_col = "SP", engine="numpyro", n.iter=20, n.burnin=10, n.thin=1)
