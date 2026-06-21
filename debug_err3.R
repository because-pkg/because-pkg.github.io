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

# First run because to get original_data and test_eqs
test_eq <- as.formula("Brain ~ Lifespan + clutch + Migration")

options(error = function() {
    sink("trace3.txt")
    traceback()
    sink()
    q("no", status=1)
})

# simulate original_data
od <- bird.dat
od$Species <- as.character(od$Species)
attr(od, "categorical_vars") <- NULL

res <- because:::run_single_dsep_test_v2(
  i = 1,
  test_eq = test_eq,
  monitor_params = "interpretable",
  engine = "jags",
  original_data = od,
  structure = bird.consensus,
  equations = eqs,
  latent = "Size",
  quiet = FALSE,
  id_col = "Species"
)
