devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
set.seed(1)
N <- 100
example_data <- data.frame(
  Y     = rpois(N, lambda = 3),
  Count = rpois(N, lambda = 5),
  Presence = rbinom(N, size = 1, prob = 0.4),
  X = rnorm(N)
)
fit_mixed.pyro <- because(
  equations = list(
    Count ~ X,
    Presence ~ X
  ),
  data = example_data,
  family = c(Count = "poisson", Presence = "binomial"),
  engine = "numpyro"
)
summary(fit_mixed.pyro)
