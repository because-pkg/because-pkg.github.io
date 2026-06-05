devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.phybase")
data(rhino.dat, package="because.phybase")
data(rhino.tree, package="because.phybase")
sem8_eq <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)

# Just run the python side directly
eq_strings <- sapply(sem8_eq, function(eq) paste(deparse(eq), collapse = ""))
for (i in seq_along(eq_strings)) {
  eq_strings[i] <- paste0(eq_strings[i], " + (1|phylo)")
}

flat_data <- flatten_for_python(rhino.dat)
flat_data$phylo <- as.integer(0:99)

py_structures <- list()
py_structures$phylo <- list(
  matrix = rhino.tree,
  transform_func = reticulate::py_run_string(numpyro_structure_definition(rhino.tree))$phylo_transform
)

py_result <- because_py$fit(
  equations = eq_strings,
  data = flat_data,
  cor_matrices = py_structures,
  num_samples = 10,
  num_warmup = 10,
  n_cores = 1
)

print(names(py_result$samples))
