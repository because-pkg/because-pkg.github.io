# Ornstein-Uhlenbeck (OU) Structural Extension for 'because'
# This script demonstrates the use of the because_structure() API.

library(because)
library(ape)

#' 1. Define the OU Kernel using the official API
#' This handles all S3 method registration automatically.
ou_kernel <- because_structure(
  name = "ou_kernel",
  precision_fn = function(tree, alpha = 0.5) {
    # Phylogenetic distance matrix
    dist_mat <- cophenetic(tree)
    n <- nrow(dist_mat)
    
    # If alpha is a vector, compute 3D array for BMA
    if (length(alpha) > 1) {
        m <- length(alpha)
        Prec_3D <- array(0, dim = c(m, n, n))
        for (j in 1:m) {
            # V_ij = exp(-alpha * d_ij)
            V <- exp(-alpha[j] * dist_mat)
            # Convert to Precision Matrix (with small ridge for stability)
            Prec_3D[j, , ] <- solve(V + diag(1e-6, n))
        }
        return(Prec_3D)
    } else {
        # Single alpha case (standard BM if alpha=0)
        V <- exp(-alpha * dist_mat)
        return(solve(V + diag(1e-6, n)))
    }
  },
  description = "Phylogenetic Ornstein-Uhlenbeck (OU) stabilizing selection kernel supporting BMA"
)

#' 2. Helper to create a BMA suite for estimating alpha
#' @param tree A phylo object.
#' @param trait_name Name of the trait for the model.
#' @param alphas Suite of selection strengths to test.
ou_structure <- function(tree, trait_name, alphas = c(0.001, 0.01, 0.1, 1.0, 10.0)) {
  # Now returns a single structure object with a 3D precision array
  # Note: trait_name is used in the model to prefix the K index (K_traitname)
  return(ou_kernel(tree = tree, alpha = alphas))
}
