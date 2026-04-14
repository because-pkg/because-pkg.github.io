devtools::load_all(".")

# Create a small dataset
set.seed(123)
n <- 100
df <- data.frame(
  X = rnorm(n),
  Y_gauss = 0.5 * rnorm(n) + 2,
  Y_poisson = rpois(n, exp(0.5 * rnorm(n)))
)

# Fit model with mixed families
# Y_gauss defaults to gaussian
# Y_poisson explicitly set to poisson
fit <- because(
  equations = list(
    Y_gauss ~ X,
    Y_poisson ~ X
  ),
  data = df,
  family = c(Y_poisson = "poisson"),
  n.iter = 1000,
  quiet = TRUE
)

# Test pp_check loop
diagnostic_nodes <- c("Y_gauss", "Y_poisson")

for (node in diagnostic_nodes) {
  cat("\nChecking node:", node, "\n")
  p <- pp_check(fit, resp = node, type = "dens_overlay")
  print(p)
  cat("Success!\n")
}
