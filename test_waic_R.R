library(because)
data(birds)

equations <- list(Brain ~ Mass)

res <- because(
    equations = equations,
    data = bird.dat,
    structure = bird.phy.cut,
    id_col = "Species",
    engine = "numpyro",
    n.iter = 50,
    n.burnin = 20,
    n.cores = 1,
    WAIC = TRUE
)

print(res$WAIC)
