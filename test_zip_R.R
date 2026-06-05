devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
set.seed(99)
N <- 200
Habitat <- rnorm(N, 0, 1)
psi <- 0.4  
lambda <- exp(0.5 + 0.7 * Habitat)
structural_zero <- rbinom(N, 1, psi)
Abundance <- ifelse(structural_zero == 1, 0, rpois(N, lambda))
zip_data <- data.frame(
  Abundance = Abundance,
  Habitat = Habitat
)

fit_zip.py <- because(
  equations = list(Abundance ~ Habitat),
  data = zip_data,
  family = c(Abundance = "zip"),
  engine="numpyro",
  n.iter=200, n.burnin=100
)

summary(fit_zip.py)
