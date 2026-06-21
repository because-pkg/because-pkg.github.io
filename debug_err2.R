devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.phybase")
library(ape)

bird.dat <- read.csv("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because_research/FinalBirds.csv", stringsAsFactors = FALSE)
bird.dat <- bird.dat[, c("Species3", "Mass", "Lifespan", "Migration", "Brain", "clutch")]
colnames(bird.dat)[1] <- "Species"
bird.dat <- na.omit(bird.dat)

bird.consensus <- read.tree("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because_research/ConsensusCut.tre")
bird.consensus <- keep.tip(bird.consensus, intersect(bird.consensus$tip.label, bird.dat$Species))
bird.dat <- bird.dat[bird.dat$Species %in% bird.consensus$tip.label, ]

eqs <- list(
  Brain ~ Size,
  Mass ~ Size,
  clutch ~ Brain,
  Migration ~ Brain,
  Lifespan ~ clutch + Migration
)

# Overwrite because to remove the tryCatch for sequential execution
# We will just redefine run_single_dsep_test_v2 to run directly
trace(because:::run_single_dsep_test_v2, tracer = quote({
  options(error = function() {
    sink("trace.txt")
    traceback()
    sink()
    q("no", status=1)
  })
}))

test_mag <- because(
  equations = eqs,
  data = bird.dat,
  structure = bird.consensus,
  id_col = "Species",
  parallel = FALSE,
  engine = "jags",
  dsep=TRUE,
  n.iter = 50,
  n.burnin = 10,
  quiet = FALSE
)
