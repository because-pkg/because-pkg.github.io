library(because)

# Simulated data
set.seed(123)
n <- 100
data <- data.frame(
    True_Mass = rnorm(n),
    julian_date = rnorm(n),
    year = factor(sample(1:5, n, replace = TRUE))
)
data$weight_g <- 10 +
    2 * data$True_Mass +
    0.5 * data$julian_date +
    rnorm(n, 0, 0.5)
data$n_pups <- rpois(n, exp(1 + 0.5 * data$True_Mass))

# Model with induced correlations (MAG)
eqs <- list(
    weight_g ~ True_Mass + julian_date,
    n_pups ~ True_Mass
)

cat("Running model and checking summary...\n")
fit <- because(
    equations = eqs,
    data = data,
    family = c(n_pups = "poisson"),
    latent = "True_Mass",
    latent_method = "correlations",
    random = list(n_pups = ~ (1 | year)),
    quiet = TRUE,
    n.iter = 1000,
    n.burnin = 500
)

# Print summary
s <- summary(fit)
print(as.data.frame(s$results))

# Check for rho
rho_found <- any(grepl("^rho", rownames(s$results)))
if (rho_found) {
    cat("\nSUCCESS: rho parameter found in summary!\n")
} else {
    cat("\nFAILURE: rho parameter NOT found in summary.\n")
    # List all params for debugging
    cat("Monitored parameters in summary:\n")
    print(rownames(s$results))
}
