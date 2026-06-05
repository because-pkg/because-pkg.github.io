library(because)
library(because.phybase)
data(rhino.dat)
data(rhino.tree)
sem8_eq <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)

fit <- because(
  equations = sem8_eq,
  data = rhino.dat,
  structure = rhino.tree,
  id_col = "SP",
  engine = "numpyro",
  n.iter=100
)
print(colnames(fit$samples[[1]]))
