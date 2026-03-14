# Test for Multiple Covariance Support (Phylo + Spatial)
# Checks if parameters are recovered correctly from a model with two additive random effects.

library(because)
library(ape)
library(MASS) # For mvrnorm
library(coda) # For summary.mcmc.list
devtools::load_all(".")

set.seed(42)

# 1. Simulation Setup
N <- 50
tree <- ape::rcoal(N)
dist_mat <- as.matrix(dist(1:N)) # Simple 1D distance for spatial
# Or strictly, make 2D spatial
coords <- matrix(runif(N * 2), N, 2)
dist_mat <- as.matrix(dist(coords))

# True Parameters
sigma_phylo <- 1.5 # tau_phylo = 1/(1.5^2) = 0.44
sigma_spatial <- 1.0 # tau_spatial = 1/(1^2) = 1.0
sigma_resid <- 0.5 # tau_resid = 1/(0.5^2) = 4.0
beta_x <- 0.8

# Covariance Matrices
V_phylo <- ape::vcv.phylo(tree)
V_spatial <- exp(-dist_mat / 0.5) # Range 0.5

# Simulate Random Effects
u_phylo <- mvrnorm(1, mu = rep(0, N), Sigma = sigma_phylo^2 * V_phylo)
u_spatial <- mvrnorm(1, mu = rep(0, N), Sigma = sigma_spatial^2 * V_spatial)
epsilon <- rnorm(N, 0, sigma_resid)

# Predictor
X <- rnorm(N)
# Response
# Y = beta*X + u_phylo + u_spatial + epsilon
Y <- beta_x * X + u_phylo + u_spatial + epsilon

data <- list(Y = Y, X = X)
equations <- list(Y ~ X)

# 2. Fit Model
# Pass list of structures
# Note: because expects Covariance-like objects (Tree or Matrix) and inverts them.
structure_list <- list(
    phylo = tree,
    spatial = V_spatial
)

message("Fitting Multi-Structure Model...")
fit <- because(
    data = data,
    equations = equations,
    structure = structure_list,
    n.iter = 5000,
    n.burnin = 1000,
    n.thin = 5,
    quiet = FALSE
)

# 3. Validation
message("\n--- Summary ---")
print(summary(fit))

# Extract posteriors from JAGS object (fit$model)
# Note: because returns 'fit' which has 'model' (rjags object) or samples?
# because returns 'mcmc.list' usually?
# Let's check because return. It returns 'result' which is mcmc.list (if parallel) or jags object?
# It returns a list with 'posteriors' (mcmc.list).

# Check variable names
params <- colnames(fit$posteriors[[1]])
message("Parameters: ", paste(params, collapse = ", "))

# Expect: tau_u_Y_phylo, tau_u_Y_spatial, tau_e_Y
# And beta_Y_X, alpha_Y

# Use summary of posteriors (mcmc.list)
summ <- summary(fit$posteriors)
get_mean <- function(param) {
    if (param %in% rownames(summ$statistics)) {
        return(summ$statistics[param, "Mean"])
    } else {
        warning(paste("Param not found:", param))
        return(NA)
    }
}

est_sigma_phylo <- get_mean("sigma_Y_phylo")
est_sigma_spatial <- get_mean("sigma_Y_spatial")
est_sigma_resid <- get_mean("sigma_Y_res")

# est_sigma_phylo is already sigma
# est_sigma_spatial is already sigma
# est_sigma_resid is already sigma

message(sprintf(
    "True Phylo Sigma: %.2f, Estimated: %.2f",
    sigma_phylo,
    est_sigma_phylo
))
message(sprintf(
    "True Spatial Sigma: %.2f, Estimated: %.2f",
    sigma_spatial,
    est_sigma_spatial
))
message(sprintf(
    "True Resid Sigma: %.2f, Estimated: %.2f",
    sigma_resid,
    est_sigma_resid
))

# Assertions
if (abs(est_sigma_phylo - sigma_phylo) > 0.5) {
    warning("Phylo sigma estimation poor")
}
if (abs(est_sigma_spatial - sigma_spatial) > 0.5) {
    warning("Spatial sigma estimation poor")
}
if (abs(est_sigma_resid - sigma_resid) > 0.5) {
    warning("Residual sigma estimation poor")
}

message("Test Complete.")
