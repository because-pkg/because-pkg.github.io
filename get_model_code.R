library(testthat); devtools::load_all(); Sys.setenv(NOT_CRAN='true'); 
context('Deterministic Logic Operators')
    set.seed(123)
    N <- 50
    A <- runif(N, 0, 10)
    in_range <- (A >= 2 & A < 5)
    Y <- A + 2 * in_range + rnorm(N, sd = 0.1)
    data <- list(Y = Y, A = A)
    eqs <- list(Y ~ A + I(A >= 2 & A < 5))
    fit <- because(equations = eqs, data = data, n.iter = 50, n.chains = 1, quiet = FALSE)
    cat('MOD_STR_START
')
    cat(fit$model_code)
    cat('
MOD_STR_END
')

