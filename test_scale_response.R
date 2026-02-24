library(devtools)
load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
set.seed(1)
df <- data.frame(
  y = scale(rnorm(100)), # 1D matrix RESPONSE
  x1 = rnorm(100)
)
eqs <- list(y ~ x1)
res <- try(because(data = df, equations = eqs))
print(res)
