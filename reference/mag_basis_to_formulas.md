# Convert basis set to because formula format

Convert basis set to because formula format

## Usage

``` r
mag_basis_to_formulas(
  basis_set,
  latent_children = NULL,
  categorical_vars = NULL,
  family = NULL,
  deterministic_terms = NULL,
  root_vars = NULL,
  hierarchical_info = NULL
)
```

## Arguments

- basis_set:

  Basis set from dagitty

- latent_children:

  Optional character vector of variables that are direct children of
  latents

- categorical_vars:

  Optional named list of categorical variable info

- family:

  Optional named character vector of family distributions

- deterministic_terms:

  Optional named list from `extract_deterministic_terms`. Internal node
  names (e.g. `"BM_x_M"`) are rendered back to their original R syntax
  (e.g. `"BM:M"`) in the output formulas, following Geiger, Verma &
  Pearl (1990).

- root_vars:

  Optional character vector of root (exogenous) node names — nodes with
  no parents in the DAG. If var1 is a root node and var2 is not, they
  are swapped so the non-root (downstream) variable becomes the
  response. This enforces the convention that parent-only variables
  should always be predictors rather than responses in independence
  tests.

## Value

List of formulas with test_var attribute
