library(testthat); devtools::load_all(); Sys.setenv(NOT_CRAN='true'); 
    set.seed(123)
    N <- 50
    A <- runif(N, 1, 10)
    Y <- 2 + 3 * A + rnorm(N, sd = 0.5)
    data <- list(Y = Y, A = A)
    eqs <- list(Y ~ A)
    fit_jags <- because(equations = eqs, data = data, engine = 'jags', n.iter = 500, n.burnin = 100, n.chains = 2, quiet = TRUE)
    fit_nimble <- because(equations = eqs, data = data, engine = 'nimble', n.iter = 500, n.burnin = 100, n.chains = 2, quiet = TRUE)
    
    cat('JAGS summary:
')
    print(fit_jags$summary$statistics)
    cat('
NIMBLE summary:
')
    print(fit_nimble$summary$statistics)

