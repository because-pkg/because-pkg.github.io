library(because)
data(storks, package="because")
equations_storks <- list(
  Storks ~ Area,
  Birth ~ Area,
  Humans ~ Birth
)
fit_storks.pyro <- because(
  equations = equations_storks,
  data = storks,
  dsep = TRUE,
  engine = "numpyro"
)
