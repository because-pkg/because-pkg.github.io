# Comprehensive Engine Parity Test: JAGS vs NIMBLE
# Verifies all supported families: Gaussian, Binomial, Poisson, Multinomial, Ordinal, ZIP, NB

library(pkgload)
load_all()
library(nimble)

# Helper for string concatenation
`%+%` <- function(a, b) paste0(a, b)

set.seed(42)
n <- 250
x <- rnorm(n)

# 1. Simulate Data for all families
# ---------------------------------
data <- data.frame(x = x)

# Gaussian
data$y_gauss <- 2 + 1.2 * x + rnorm(n, 0, 1)

# Binomial
data$y_bin <- rbinom(n, 1, plogis(-0.5 + 0.8 * x))

# Poisson
data$y_pois <- rpois(n, exp(0.5 + 0.5 * x))

# Negative Binomial (NB)
# mu = exp(0.5 + 0.5*x), size = 2
data$y_nb <- rnbinom(n, size = 2, mu = exp(0.5 + 0.5 * x))

# Zero-Inflated Poisson (ZIP)
# psi = 0.3 (prob of zero), mu = exp(1 + 0.5*x)
is_zero <- rbinom(n, 1, 0.3)
data$y_zip <- ifelse(is_zero == 1, 0, rpois(n, exp(1 + 0.5 * x)))

# Multinomial (3 categories)
# eta2 = -0.5 + 1*x, eta3 = 0.5 - 0.5*x, eta1 = 0
eta2 <- -0.5 + 1 * x
eta3 <- 0.5 - 0.5 * x
probs <- t(apply(cbind(0, eta2, eta3), 1, function(e) exp(e) / sum(exp(e))))
data$y_multinom <- apply(probs, 1, function(p) sample(1:3, 1, prob = p))

# Ordinal (3 categories)
# eta = 0.8 * x, cutpoints = -1, 1
eta_ord <- 0.8 * x
p1 <- plogis(-1 - eta_ord)
p2 <- plogis(1 - eta_ord) - p1
p3 <- 1 - plogis(1 - eta_ord)
probs_ord <- cbind(p1, p2, p3)
data$y_ordinal <- apply(probs_ord, 1, function(p) sample(1:3, 1, prob = p))

# 2. Define Test Scenarios
# ------------------------
models <- list(
    gaussian = list(eqs = list(y_gauss ~ x), family = c(y_gauss = "gaussian")),
    binomial = list(eqs = list(y_bin ~ x), family = c(y_bin = "binomial")),
    poisson = list(eqs = list(y_pois ~ x), family = c(y_pois = "poisson")),
    negbin = list(eqs = list(y_nb ~ x), family = c(y_nb = "negbinomial")),
    zip = list(eqs = list(y_zip ~ x), family = c(y_zip = "zip")),
    multinom = list(
        eqs = list(y_multinom ~ x),
        family = c(y_multinom = "multinomial")
    ),
    ordinal = list(eqs = list(y_ordinal ~ x), family = c(y_ordinal = "ordinal"))
)

results <- list()

# 3. Fit and Compare
# ------------------
for (name in names(models)) {
    cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
    cat(sprintf("Testing Family: %s\n", toupper(name)))
    cat(paste(rep("=", 60), collapse = "") %+% "\n")

    mod <- models[[name]]

    cat("Fitting JAGS...\n")
    fit_jags <- because(
        mod$eqs,
        data = data,
        family = mod$family,
        engine = "jags",
        n.iter = 5000,
        n.burnin = 1000,
        quiet = TRUE
    )

    cat("Fitting NIMBLE...\n")
    fit_nimble <- because(
        mod$eqs,
        data = data,
        family = mod$family,
        engine = "nimble",
        n.iter = 5000,
        n.burnin = 1000,
        quiet = TRUE
    )

    jags_sum <- as.data.frame(fit_jags$summary$statistics)
    nimble_sum <- as.data.frame(fit_nimble$summary$statistics)

    params <- intersect(rownames(jags_sum), rownames(nimble_sum))
    # Standardize parameter selection: include intercepts, coefficients, cutpoints, 
    # correlation terms, and dispersion/inflation parameters
    params <- params[grepl("beta_|alpha_|cutpoint_|rho_|r_|psi_", params)]

    comp <- data.frame(
        Param = params,
        JAGS_Mean = round(jags_sum[params, "Mean"], 3),
        NIMBLE_Mean = round(nimble_sum[params, "Mean"], 3),
        JAGS_SD = round(jags_sum[params, "SD"], 3),
        NIMBLE_SD = round(nimble_sum[params, "SD"], 3)
    )

    # Parity Metric: Difference in units of Posterior SD
    comp$Diff_SD <- round(
        abs(comp$JAGS_Mean - comp$NIMBLE_Mean) / comp$JAGS_SD,
        3
    )
    comp$Status <- ifelse(comp$Diff_SD < 0.6, "PASS", "CHECK")

    print(comp[, c("Param", "JAGS_Mean", "NIMBLE_Mean", "Diff_SD", "Status")])
    results[[name]] <- comp
}

cat("\n" %+% paste(rep("-", 40), collapse = "") %+% "\n")
cat("SUMMARY OF PARITY CHECKS\n")
cat(paste(rep("-", 40), collapse = "") %+% "\n")
for (name in names(results)) {
    status <- if (all(results[[name]]$Status == "PASS")) {
        "PASSED"
    } else {
        "NEEDS REVIEW"
    }
    cat(sprintf(" - %-10s: %s\n", name, status))
}
