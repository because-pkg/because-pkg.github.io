devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because.phybase")
py_def <- numpyro_structure_definition.phylo()
env <- reticulate::py_run_string(py_def)
py_structures <- list()
py_structures[["phylo"]] <- list(
  matrix = matrix(1:4, 2, 2),
  transform_func = env[["phylo_transform"]]
)

py_tester <- reticulate::py_run_string("
def check_type(obj):
    print(f'Type of obj is: {type(obj)}')
    if isinstance(obj, dict):
        print(f'Keys: {obj.keys()}')
        if 'phylo' in obj:
            print(f'Type of phylo: {type(obj[\"phylo\"])}')
            if isinstance(obj[\"phylo\"], dict):
                print(f'Keys of phylo: {obj[\"phylo\"].keys()}')
")

py_tester$check_type(py_structures)
