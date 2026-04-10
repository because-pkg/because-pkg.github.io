## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## -----------------------------------------------------------------------------
# # load data
# library(because)
# data(storks)
# plot(storks$Storks, storks$Birth,
#   xlab = "Number of Stork Pairs",
#   ylab = "Number of Babies Born (thousands)"
# )

## ----echo=FALSE, eval=TRUE, out.width="100%"----------------------------------
knitr::include_graphics("figures/storks_vs_babies.png")

## -----------------------------------------------------------------------------
# lm_storks <- lm(Birth ~ Storks, data = storks)
# summary(lm_storks)

## -----------------------------------------------------------------------------
# birth_area.lm <- lm(Birth ~ Area, data = storks)
# summary(birth_area.lm)
# # Coefficients:
# #              Estimate Std. Error t value Pr(>|t|)
# # (Intercept) -7.7754992 56.9376784  -0.137    0.893
# # Area         0.0017229  0.0001861   9.259 1.36e-07 ***
# 
# storks_area.lm <- lm(Storks ~ Area, data = storks)
# summary(storks_area.lm)
# # Coefficients:
# #              Estimate Std. Error t value Pr(>|t|)
# # (Intercept) -6.069e+01  2.591e+03  -0.023   0.9816
# # Area         2.331e-02  8.467e-03   2.753   0.0148 *

## -----------------------------------------------------------------------------
# lm_birth_storks_area <- lm(Birth ~ Storks + Area,
#   data = storks
# )
# summary(lm_birth_storks_area)
# 
# # Coefficients:
# #               Estimate Std. Error t value Pr(>|t|)
# # (Intercept) -7.4116870 56.7021798  -0.131    0.898
# # Storks       0.0059949  0.0056510   1.061    0.307
# # Area         0.0015832  0.0002273   6.964 6.62e-06 ***

## -----------------------------------------------------------------------------
# library(because)
# dag_storks <- list(
#   Storks ~ Area,
#   Birth ~ Area)
# 
# plot_dag(dag_storks)

## ----echo=FALSE, eval=TRUE, out.width="100%"----------------------------------
knitr::include_graphics("figures/dag_storks.png")

## -----------------------------------------------------------------------------
# dag_storks2 <- list(Storks ~ Area,
#                     Birth ~ Area,
#                     Humans ~ Birth)
# 
# plot_dag(dag_storks2)
# 

## ----echo=FALSE, eval=TRUE, out.width="100%"----------------------------------
knitr::include_graphics("figures/dag_storks2.png")

## -----------------------------------------------------------------------------
# equations_storks <- list(
#   Storks ~ Area,
#   Birth ~ Area,
#   Humans ~ Birth
# )

## -----------------------------------------------------------------------------
# storks <- scale(storks[2:5])

## -----------------------------------------------------------------------------
# fit_storks <- because(
#   equations = equations_storks,
#   data = storks,
#   dsep = TRUE
# )
# summary(fit_storks)

## -----------------------------------------------------------------------------
# fit_storks_final <- because(
#   equations = equations_storks,
#   data = storks
# )
# summary(fit_storks_final)

## -----------------------------------------------------------------------------
# plot_dag(fit_storks_final)

