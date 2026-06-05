devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
set.seed(42)
N <- 200

# Predictors
Temperature <- rnorm(N, mean = 15, sd = 5)
Precipitation <- rnorm(N, mean = 800, sd = 200)

# True logistic relationship
logit_p <- -1 + 0.15 * scale(Temperature) + 0.08 * scale(Precipitation)
Presence <- rbinom(N, size = 1, prob = 1 / (1 + exp(-logit_p)))

presence_data <- data.frame(
  Presence = Presence,
  Temperature = scale(Temperature)[, 1],
  Precipitation = scale(Precipitation)[, 1]
)
fit_bern <- because(
  equations = list(Presence ~ Temperature + Precipitation),
  data = presence_data,
  family = c(Presence = "bernoulli"),
  n.iter = 100, n.burnin = 50
)

summary(fit_bern)
