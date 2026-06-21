library(because)
dag_adj <- matrix(0, nrow=6, ncol=6)
rownames(dag_adj) <- colnames(dag_adj) <- c("Size", "Brain", "Mass", "Migration", "Clutch", "Lifespan")
dag_adj["Size", "Brain"] <- 1
dag_adj["Size", "Mass"] <- 1
dag_adj["Brain", "Clutch"] <- 1
dag_adj["Brain", "Migration"] <- 1
dag_adj["Clutch", "Lifespan"] <- 1
dag_adj["Migration", "Lifespan"] <- 1

basis <- because:::because_dsep(dag_adj, latent="Size")
for (f in basis$tests) {
  print(f)
}
