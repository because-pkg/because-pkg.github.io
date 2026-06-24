library(because)
set.seed(123)
data <- data.frame(A = rnorm(100), B = rnorm(100), C = rnorm(100))
fit <- because(list(B ~ A, C ~ B), data = data, n.iter = 1000, dsep = TRUE, silent = TRUE)

s95 <- summary(fit, rope = c(-0.1, 0.1), eq_test = TRUE, prob = 0.95)
print(s95$results)

s99 <- summary(fit, rope = c(-0.1, 0.1), eq_test = TRUE, prob = 0.99)
print(s99$results)
