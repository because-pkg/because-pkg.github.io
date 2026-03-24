# Format Data for Because Analysis

Converts data from long format (one row per observation) to the list
format required by
[`because`](https://because-pkg.github.io/because/reference/because.md).

## Usage

``` r
because_format_data(data, species_col = "SP", tree = NULL)
```

## Arguments

- data:

  A data.frame in long format with one row per observation.

- species_col:

  Name of the column containing species or unit identifiers (default:
  "SP").

- tree:

  A phylogenetic tree or other structural object. Optional. If provided,
  determines the order of units.

## Value

A named list where each element is either:

- A numeric vector (if all species have exactly 1 observation)

- A numeric matrix with species in rows and replicates in columns

Species are ordered to match `tree$tip.label` (if provided) or sorted
alphabetically.

## Details

This function handles:

- Different numbers of replicates per species (creates rectangular
  matrix with NA padding)

- Missing values (NA)

- Automatic alignment with phylogenetic tree tip labels (if provided)

When species have different numbers of replicates, the function creates
a matrix with dimensions (number of species) x (maximum number of
replicates). Species with fewer replicates are padded with NA values.

If a tree is provided: Species in the tree but not in the data will have
all NA values. Species in the data but not in the tree will be excluded
with a warning.

If no tree is provided: All species in the data are included, sorted
alphabetically by their ID.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example data in long format
data_long <- data.frame(
  SP = c("sp1", "sp1", "sp1", "sp2", "sp2", "sp3"),
  BM = c(1.2, 1.3, 1.1, 2.1, 2.2, 1.8),
  NL = c(0.5, 0.6, NA, 0.7, 0.8, 0.9)
)

# With tree
if (requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::read.tree(text = "(sp1:1,sp2:1,sp3:1);")
  data_list <- because_format_data(data_long, species_col = "SP", tree = tree)
}

# Without tree (general repeated measures)
data_list_no_tree <- because_format_data(data_long, species_col = "SP")
} # }
```
