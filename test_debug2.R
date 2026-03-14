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

cat("\nExtracting generated JAGS model via simple because_model call:\n")
prep_data <- because:::prepare_hierarchical_jags_data(
  list(individual = dataset, enviro = density),
  hierarchical_info$levels,
  hierarchical_info$hierarchy,
  hierarchical_info$link_vars
)
hier_info_for_jags <- list(
  levels = hierarchical_info$levels,
  hierarchy = hierarchical_info$hierarchy,
  link_vars = hierarchical_info$link_vars,
  data = list(individual = dataset, enviro = density)
)

mod <- because:::because_model(
  equations = list(julian_date ~ mean_adsub_size + BGI + EOS),
  random_terms = random_terms,
  hierarchical_info = hier_info_for_jags,
  data_list = prep_data$jags_data
)

cat(paste(mod$model_lines, collapse = "\n"))
