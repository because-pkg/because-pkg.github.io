devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.phybase")
library(ape)

eqs <- list(Brain ~ Lifespan + clutch + Migration)
latent <- "Size"
random_terms <- list()

trace(because:::because_dsep, tracer = quote({
  if (!dsep) {
    print("INSIDE because_dsep")
    # We will let it run, and print basis before mag_basis_to_formulas
  }
}), at = 1) # wait, line number trace is tricky

# let's just override because_dsep locally
# actually, let's just run because_dsep
res <- because:::because_dsep(
  equations = eqs,
  latent = latent,
  random_terms = random_terms,
  quiet = FALSE
)
print(res)
