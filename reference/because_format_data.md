# Format Long-Format Repeated Measures Data for Because Analysis

Converts data from long format (multiple rows per unit/group, such as
repeated measurements per individual, site survey visits, or
within-species replicates) to the named list format (unit-level vectors
and replicate matrices) used by
[`because`](https://because-pkg.github.io/because/reference/because.md)
for error-in-variables (`variability = "reps"`), detection/occupancy
models, and structured covariance models.

## Usage

``` r
because_format_data(
  data,
  id_col = "ID",
  structure = NULL,
  species_col = NULL,
  tree = NULL
)
```

## Arguments

- data:

  A `data.frame` in long format with one row per observation.

- id_col:

  Name of the column containing unit identifiers (e.g., individual,
  site, or species IDs; default: `"ID"`, falling back to `"SP"` if
  present).

- structure:

  Optional covariance structure object (e.g., a spatial/pedigree matrix
  with row names, a `because_structure` object, or a phylogenetic tree).
  If provided, determines the ordering of units in the output.

- species_col:

  Deprecated alias for `id_col` (retained for backward compatibility).

- tree:

  Deprecated alias for `structure` (retained for backward
  compatibility).

## Value

A named list where each element is either:

- A named numeric vector (if a variable is constant within each unit or
  all units have 1 observation)

- A numeric matrix with units in rows and replicates in columns (padded
  with `NA`)

Units are ordered to match the labels of `structure` (if provided) or
sorted alphabetically.

## Details

This function handles:

- Different numbers of replicates per unit (creates a rectangular matrix
  with `NA` padding)

- Automatic distinction between unit-constant covariates (returned as
  vectors) and replicate-varying measurements (returned as matrices)

- Missing values (`NA`)

- Automatic alignment with covariance structure labels (e.g., matrix row
  names or phylogenetic tree tip labels, if `structure` is provided)

When units have different numbers of replicates, the function creates a
matrix with dimensions
`(number of units) x (maximum number of replicates)`. Units with fewer
replicates are padded with `NA` values.

If a `structure` is provided:

- Units present in the structure but missing from `data` receive all
  `NA` values (enabling Bayesian imputation).

- Units present in `data` but absent from the structure are excluded
  with a warning.

If no `structure` is provided: All units in `data` are included, sorted
alphabetically by their identifier.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: General repeated measures / Error-in-Variables (no structure)
data_long <- data.frame(
  ID   = c("ind1", "ind1", "ind1", "ind2", "ind2", "ind3"),
  Age  = c(3, 3, 3, 5, 5, 2),          # Constant within individual -> vector
  Mass = c(1.2, 1.3, 1.1, 2.1, 2.2, 1.8), # Varying -> matrix (3 x 3)
  Cort = c(0.5, 0.6, NA, 0.7, 0.8, 0.9)   # Varying -> matrix (3 x 3)
)
data_list <- because_format_data(data_long, id_col = "ID")

# Example 2: Aligning with a covariance matrix (e.g., spatial or pedigree)
prec_mat <- diag(3)
rownames(prec_mat) <- colnames(prec_mat) <- c("ind3", "ind1", "ind2")
data_aligned <- because_format_data(data_long, id_col = "ID", structure = prec_mat)
} # }
```
