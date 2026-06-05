library(reticulate)
because_py <- import("because")
family <- c(Count = "poisson", Presence = "binomial")
print(str(r_to_py(family)))
