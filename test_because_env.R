devtools::load_all()

cat("\nTesting because() auto-binding...\n")
set.seed(42)
N <- 100
x <- rnorm(N, 0, 1)
y_counts <- rpois(N, exp(0.5 + 0.8 * x))
z <- rnorm(N, 0.2 + 0.4 * x, 1)

data_dict <- data.frame(x = x, z = z, y = y_counts)

# Provide formulas as a list
fit <- because(
  equations = list(y ~ x, z ~ x),
  data = data_dict,
  family = list(y = "poisson"),
  engine = "numpyro",
  n.iter = 500,
  n.chains = 1
)

print(summary(fit))
cat("\nSUCCESS! The package automatically bound to because_env and ran the model.\n")
