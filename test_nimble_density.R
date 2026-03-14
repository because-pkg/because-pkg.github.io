library(nimble)

model_string <- "{
  for(i in 1:10) {
    # Testing dnorm in expression
    ld_norm[i] <- dnorm(y[i], mu, 1, log = 1)
    # Testing dpois in expression
    ld_pois[i] <- dpois(counts[i], lambda, log = 1)
    
    # Example likelihood trick
    total_lik[i] <- exp(ld_norm[i]) + exp(ld_pois[i])
    
    # dummy stochastic
    z[i] ~ dbern(0.5)
  }
  mu ~ dnorm(0, 1)
  lambda ~ dgamma(1, 1)
}"

code <- parse(text = model_string)[[1]]

tryCatch({
  m1 <- nimbleModel(code = code, constants=list(y=rnorm(10), counts=rpois(10, 2)))
  print("nimbleModel succeeded!")
}, error = function(e) print(paste("nimbleModel failed:", e$message)))

