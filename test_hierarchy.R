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
  name = paste0("N", 1:(n_enviro * n_indiv_per_enviro)),
  burrow_system = sample(1:5, n_enviro * n_indiv_per_enviro, replace = TRUE)
)

w_dag_full7_a <- list(
  weight_g ~ EOS + ageclass + BGI + julian_date + mean_adsub_size + sex,
  mean_adsub_size ~ BGI,
  julian_date ~ BGI + EOS,
  EOS ~ L1,
  BGI ~ L1
)

# Run dsep directly
dsep_res <- because_dsep(
  equations = w_dag_full7_a,
  latent = c("L1"),
  hierarchical_info = list(
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
    hierarchy = "enviro > individual"
  ),
  quiet = FALSE
)

# Actually run the model
cat("\nRunning because model with hierarchy dsep true...\n")
fit <- because(
  data = list(individual = dataset, enviro = density),
  equations = w_dag_full7_a,
  latent = c("L1"),
  random = ~ (1 | name) + (1 | burrow_system),
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
  link_vars = c("year", "zone"),
  dsep = TRUE,
  n.iter = 100,
  n.burnin = 10,
  parallel = FALSE,
  quiet = FALSE
)

cat("Successfully ran dsep tests!\n")
