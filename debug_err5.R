devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.phybase")
library(ape)

eqs <- list(Brain ~ Lifespan + clutch + Migration)
latent <- "Size"
random_terms <- list()

options(error = function() {
    sink("trace5.txt")
    traceback()
    sink()
    q("no", status=1)
})

trace(because:::mag_basis_to_formulas, tracer = quote({
  print(str(basis_set))
}))

res <- because:::because_dsep(
  equations = eqs,
  latent = latent,
  random_terms = random_terms,
  quiet = FALSE
)
