## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>", eval = FALSE
)

## -----------------------------------------------------------------------------
# library(because)
# 
# # Define your SEM
# equations <- list(
#     Mass ~ Age + Sex,
#     Lifespan ~ Mass + Habitat
# )
# 
# # Fit using NIMBLE instead of JAGS
# fit <- because(
#     equations = equations,
#     data = my_data,
#     engine = "nimble"
# )

## -----------------------------------------------------------------------------
# fit <- because(
#     equations = equations,
#     data = my_data,
#     engine = "nimble", # Works for "jags" too!
#     parallel = TRUE,
#     n.chains = 4,
#     n.cores = 4
# )

## -----------------------------------------------------------------------------
# # Override the default RW-MH sampler for a specific coefficient
# fit <- because(
#     equations = equations,
#     data = my_data,
#     engine = "nimble",
#     nimble_samplers = list(
#         "beta_Lifespan_Mass" = "slice"
#     )
# )

