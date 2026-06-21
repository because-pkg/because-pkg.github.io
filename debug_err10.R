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

# Use trace to print dsep_equations_to_run right before calling python
trace(because:::because, tracer = quote({
  if (engine == "numpyro" && dsep) {
    cat("Num dsep equations to run from R:", length(dsep_equations_to_run), "\n")
  }
}), at = 2406)

tryCatch({
test_mag <- because(
  equations = eqs,
  data = bird.dat,
  structure = bird.consensus,
  id_col = "Species",
  parallel = FALSE,
  engine = "numpyro",
  dsep = TRUE,
  n.iter = 5,
  n.burnin = 1,
  quiet = FALSE
)
}, error=function(e) print(e))
