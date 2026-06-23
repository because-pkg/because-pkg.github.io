library(coda)
m <- as.mcmc(matrix(rnorm(1000), ncol=2))
s <- summary(m, quantiles=c(0.005, 0.5, 0.995))
print(colnames(s$quantiles))
