library(because)
library(because.phybase)
library(ape)

data(rhino.dat)
data(rhino.tree)

# Create a multiPhylo object by adding some random noise to branch lengths
multi_tree <- rmtree(10, 100) # 10 trees, 100 tips
# Give them the same tip labels as rhino.tree
for(i in 1:10) {
  multi_tree[[i]]$tip.label <- rhino.tree$tip.label
}

sem8_eq <- list(
    LS ~ BM,
    NL ~ BM + RS,
    DD ~ NL
)

cat("\n=== RUNNING NUMPYRO MULTIPHYLO ===\n")
t1 <- Sys.time()
fit_numpyro <- because(
    equations = sem8_eq,
    data = rhino.dat,
    structure = multi_tree,
    id_col = "SP",
    engine = "numpyro",
    parallel = TRUE,
    n.cores = 3,
    n.iter = 1000,
    n.warmup = 500
)
t2 <- Sys.time()
cat(sprintf("NumPyro time: %.2f secs\n", as.numeric(difftime(t2, t1, units="secs"))))


cat("\n=== RUNNING NIMBLE MULTIPHYLO ===\n")
t1 <- Sys.time()
fit_nimble <- because(
    equations = sem8_eq,
    data = rhino.dat,
    structure = multi_tree,
    id_col = "SP",
    engine = "nimble",
    parallel = TRUE,
    n.cores = 3,
    n.iter = 1000,
    n.burnin = 500
)
t2 <- Sys.time()
cat(sprintf("NIMBLE time: %.2f secs\n", as.numeric(difftime(t2, t1, units="secs"))))
