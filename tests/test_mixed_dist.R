library(because)

# Reproduce mixed distribution bug
set.seed(123)
N <- 40
df <- data.frame(
    SP = paste0("sp", 1:N),
    Y_cont = rnorm(N, mean = 10, sd = 2),
    Y_count = rpois(N, lambda = 5),
    X = rnorm(N)
)

cat("Testing mixed distribution model...\n")
fit <- because(
    data = df,
    id_col = "SP",
    equations = list(
        Y_cont ~ X,
        Y_count ~ X
    ),
    family = list(
        Y_cont = "gaussian",
        Y_count = "poisson"
    ),
    WAIC = TRUE,
    n.chains = 2,
    n.iter = 500,
    quiet = FALSE
)

cat("\n✓ SUCCESS!\n")
print(fit$WAIC)
