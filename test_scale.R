library(because)
set.seed(1)
df <- data.frame(
  y = rnorm(100),
  x1 = as.numeric(scale(rnorm(100))), # proper vector
  x2 = scale(rnorm(100)) # 1D matrix
)
eqs <- list(y ~ x1 + x2)
# should fail or succeed?
res <- because(data = df, equations = eqs)
print(res)
