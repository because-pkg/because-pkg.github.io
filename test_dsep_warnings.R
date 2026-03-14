library(because)
# Source all R files to ensure local changes are used
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
    source(f)
}

# Create minimal hierarchical data with categorical variables
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
dataset$ageclass <- as.factor(dataset$ageclass)
# dataset$sex <- as.factor(dataset$sex) # Keep as numeric 0/1 for binomial

# Smaller DAG for faster testing
# Test 1: continuous ~ categorical (julian_date ~ ageclass)
# Test 2: categorical ~ continuous (ageclass ~ sex) -> Wait, sex also factor
# Test 3: categorical ~ continuous (ageclass ~ EOS)

w_dag <- list(
    weight_g ~ EOS + ageclass + BGI + julian_date + sex,
    julian_date ~ BGI + EOS,
    EOS ~ L1,
    BGI ~ L1
)

cat("\nRunning because with dsep=TRUE...\n")
model <- because(
    data = list(
        individual = dataset,
        enviro = density
    ),
    equations = w_dag,
    random = ~ (1 | name) + (1 | burrow_system),
    latent = c("L1"),
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
    family = c(ageclass = "multinomial", sex = "binomial"),
    dsep = TRUE,
    n.iter = 10,
    n.adapt = 10,
    n.burnin = 0, # Very fast for just testing parameter existence
    parallel = FALSE
)

cat("\nCalling summary(model)...\n")
# cat("\nGenerated JAGS Model:\n", model$model, "\n")
s <- summary(model)
print(s)
