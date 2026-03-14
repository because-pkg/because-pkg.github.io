library(because)
library(ape)
devtools::load_all()

set.seed(123)
n_species <- 30
tree <- rtree(n_species)
tree$tip.label <- paste0("sp", 1:n_species)

# Simulate phylogenetic signal
# We use a simple Brownian motion process
vcv_mat <- vcv(tree)
vcv_mat <- vcv_mat / max(vcv_mat) # standardize

# Cholesky decomposition for simulation
L <- t(chol(vcv_mat))
u <- L %*% rnorm(n_species)

# Predictor
X <- rnorm(n_species)
# Response with phylogenetic signal + predictor effect + noise
Y <- 2 + 1.5 * X + as.numeric(u) * 2 + rnorm(n_species, 0, 0.5)

dat <- data.frame(
  species = tree$tip.label,
  Y = Y,
  X = X
)

eqs <- list(Y ~ X)

message("Fitting PHYLOGENETIC model with NIMBLE...")
fit_nimble_phylo <- because(
  equations = eqs,
  data = dat,
  tree = tree,
  engine = 'nimble',
  n.iter = 100,
  n.burnin = 10,
  n.chains = 2,
  quiet = FALSE
)

print(fit_nimble_phylo$summary$statistics)
message("Done!")
