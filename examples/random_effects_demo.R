# because Random Effects Demonstration
# This script demonstrates how to use the new lmer-style random effects syntax
# (e.g., (1|Group)) in because models.

library(because)

set.seed(42)

# --- Example 1: Gaussian Random Intercept Model ---
message("\n--- Example 1: Gaussian Random Intercept ---")

# Simulate Data
N <- 100
n_groups <- 10
groups <- sample(LETTERS[1:n_groups], N, replace = TRUE)
group_eff <- rnorm(n_groups, 0, 1.5)
names(group_eff) <- LETTERS[1:n_groups]

x <- rnorm(N)
y <- 1 + 0.5 * x + group_eff[groups] + rnorm(N, 0, 1)
df <- data.frame(Y = y, X = x, Group = as.factor(groups))

# Run because with random effects syntax
# Note: optimise=TRUE is automatic when random effects are detected without tree
fit_gaussian <- because(
    equations = list(Y ~ X + (1 | Group)),
    data = df,
    n.chains = 1, # Short run for demo
    n.iter = 1000,
    n.burnin = 200,
    quiet = TRUE
)

# Inspect Results
print(fit_gaussian$summary$statistics[
    c("alphaY", "beta_Y_X", "sigma_Y_Group", "sigma_Y_res"),
    "Mean"
])
message("Expected: alphaY~1, beta~0.5, sigma_Y_Group~1.5, sigma_Y_res~1")


# --- Example 2: Poisson GLMM with Random Intercept ---
message("\n--- Example 2: Poisson GLMM ---")

# Simulate Poisson Data with Overdispersion/random effects
lambda_log <- 0.5 + 0.3 * x + group_eff[groups]
y_pois <- rpois(N, exp(lambda_log))
df_pois <- data.frame(Y = y_pois, X = x, Group = as.factor(groups))

fit_poisson <- because(
    equations = list(Y ~ X + (1 | Group)),
    data = df_pois,
    distribution = list(Y = "poisson"),
    n.chains = 1,
    n.iter = 1000,
    n.burnin = 200,
    quiet = TRUE
)

# Inspect Results
# Note: For Poisson, sigma_Y_res corresponds to overdispersion (epsilon)
print(fit_poisson$summary$statistics[
    c("alphaY", "beta_Y_X", "sigma_Y_Group"),
    "Mean"
])


# --- Example 3: Multiple Random Effects ---
message("\n--- Example 3: Multiple Random Effects ---")
# Y ~ X + (1|Site) + (1|Year)

n_sites <- 5
n_years <- 3
sites <- sample(paste0("S", 1:n_sites), N, replace = TRUE)
years <- sample(paste0("Y", 1:n_years), N, replace = TRUE)

df_multi <- data.frame(Y = y, X = x, Site = sites, Year = years)

fit_multi <- because(
    equations = list(Y ~ X + (1 | Site) + (1 | Year)),
    data = df_multi,
    n.chains = 1,
    n.iter = 500, # Very short run
    quiet = TRUE
)

print(rownames(fit_multi$summary$statistics))
# Should see sigma_Y_Site and sigma_Y_Year
