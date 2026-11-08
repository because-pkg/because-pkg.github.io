# Process covariance structures for because models

Validates and formats user-supplied covariance structures (spatial,
phylogenetic, etc.) into the data and metadata needed by
JAGS/NIMBLE/NumPyro.

## Usage

``` r
process_because_structures(
  structure,
  structures,
  structure_obj,
  data,
  hierarchical_info,
  is_hierarchical,
  equations,
  family,
  family_obj,
  levels,
  multiscale,
  latent,
  quiet,
  engine = "jags",
  original_raw_data_before_preprocess = NULL,
  row_ids = NULL
)
```
