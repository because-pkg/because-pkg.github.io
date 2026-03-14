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

dsep_res <- because:::because_dsep(
  equations = list(
    weight_g ~ EOS + ageclass + BGI + julian_date + mean_adsub_size + sex,
    mean_adsub_size ~ BGI,
    julian_date ~ BGI + EOS,
    EOS ~ L1,
    BGI ~ L1
  ),
  latent = c("L1"),
  hierarchical_info = list(
    data = list(individual = dataset, enviro = density),
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
  ),
  quiet = FALSE
)

cat("Tests generated:\n")

# Try because_model on the first test
test_eq <- dsep_res$tests[[1]]
random_terms <- because:::parse_global_random(
  ~ (1 | name) + (1 | burrow_system),
  list(test_eq)
)

cat(
  "\nExtracting generated JAGS model via simple because_model call for test 1:\n"
)

hierarchical_info <- list(
  data = list(individual = dataset, enviro = density),
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

prep_res <- because:::prepare_hierarchical_jags_data(
  hierarchical_info,
  c("julian_date", "mean_adsub_size", "BGI", "EOS", "name", "burrow_system")
)
hierarchical_info$data <- prep_res$data_for_r

mod <- because:::because_model(
  equations = list(test_eq),
  random_terms = random_terms,
  hierarchical_info = hierarchical_info
)
cat("\nJAGS MODEL LINES:\n")
cat(paste(mod$model_lines, collapse = "\n"))
