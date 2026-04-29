# Auto-Stack Multispecies Data

Converts a list of matrices (Wide) into a stacked dataframe (Long) for
multiscale analysis. Replicates site and species covariates accordingly.

## Usage

``` r
auto_stack_multispecies_data(data, equations, quiet = FALSE)
```

## Arguments

- data:

  Input list of data (e.g. list(Y=list(Sp1=mat...), Hab=vec, Trait=vec))

- equations:

  List of model equations

- quiet:

  Logical, suppress messages

## Value

List with components:

- data:

  New stacked dataframe

- random_part:

  String to append to random formula (e.g. "+ (1\|SpeciesID)")

- is_stacked:

  Logical, whether stacking occurred
