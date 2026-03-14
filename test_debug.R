pkgload::load_all()

# Create minimal hierarchical data
set.seed(123)
n_enviro <- 5
n_indiv_per_enviro <- 10

density <- data.frame(
  zone = 1:n_enviro,
  year = rep(2020, n_enviro),
  EOS = rnorm(n_enviro),
  BGI = rnorm(n_enviro),
  mean_adsub_size = rnorm(n_enviro)
)

dataset <- data.frame(
  id = 1:(n_enviro * n_indiv_per_enviro),
  zone = rep(1:n_enviro, each = n_indiv_per_enviro),
  year = rep(2020, n_enviro * n_indiv_per_enviro),
  ageclass = sample(0:3, n_enviro * n_indiv_per_enviro, replace = TRUE),
  sex = sample(0:1, n_enviro * n_indiv_per_enviro, replace = TRUE),
  weight_g = rnorm(n_enviro * n_indiv_per_enviro),
  julian_date = runif(n_enviro * n_indiv_per_enviro),
  name = paste0("name", 1:(n_enviro * n_indiv_per_enviro)),
  burrow_system = sample(1:5, n_enviro * n_indiv_per_enviro, replace = TRUE)
)

data_list = list(individual = dataset, enviro = density)
hierarchical_info <- list(
  data = data_list,
  levels = list(
    individual = c(
      "ageclass",
      "burrow_system",
      "julian_date",
      "weight_g",
      "name",
      "sex",
      "year"
    ),
    enviro = c("BGI", "EOS", "mean_adsub_size")
  ),
  hierarchy = "enviro > individual",
  link_vars = c("year", "zone")
)
# Parse the random effects for julian_date
random_terms <- because:::parse_global_random(
  ~ (1 | name) + (1 | burrow_system),
  list(julian_date ~ mean_adsub_size + BGI + EOS)
)
print("Random Terms parsed:")
print(random_terms)

# Generate JAGS output to see what the loop bound for burrow_system is
cat("\nTesting model script generation directly\n")
fit <- because(
  data = list(individual = dataset, enviro = density), # need original nested list format
  equations = list(julian_date ~ mean_adsub_size + BGI + EOS),
  random = ~ (1 | name) + (1 | burrow_system),
  levels = hierarchical_info$levels,
  hierarchy = hierarchical_info$hierarchy,
  link_vars = hierarchical_info$link_vars,
  n.iter = 1,
  n.adapt = 1,
  n.burnin = 0 # just compile/run minimal
)
