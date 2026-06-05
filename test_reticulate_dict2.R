library(reticulate)
family <- c(Count = "poisson", Presence = "binomial")
print(str(r_to_py(as.list(family))))
