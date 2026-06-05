devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.phybase")
data(rhino.dat, package="because.phybase")
data(rhino.tree, package="because.phybase")
sem8_eq <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)

fit <- because(equations = sem8_eq, data = rhino.dat, structure = rhino.tree, id_col = "SP", engine = "jags", n.iter=100)
print(fit)
