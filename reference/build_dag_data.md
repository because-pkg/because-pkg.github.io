# Build tidy DAG data frame for a single model or formula list

Build tidy DAG data frame for a single model or formula list

## Usage

``` r
build_dag_data(
  obj,
  label,
  latent = NULL,
  family = NULL,
  layout = "kk",
  coords = NULL,
  node_size = 12,
  text_size = 3.5,
  edge_label_size = 3,
  edge_color_scheme = "directional",
  type = "raw",
  multinomial_probabilities = TRUE,
  show_coefficients = TRUE
)
```

## Arguments

- obj:

  A `because` model object or a list of formulas.

- label:

  Character label for this model (used in faceting).

- latent:

  Optional character vector of latent variable names.

- family:

  Optional named character vector of family distributions.

- layout:

  Layout algorithm (passed to ggdag::tidy_dagitty).

- coords:

  Optional named list of node coordinates.

- node_size:

  Base node size.

- text_size:

  Node label text size.

- edge_label_size:

  Edge coefficient label size.

- edge_color_scheme:

  One of "directional", "binary", "monochrome".

- type:

  "raw" or "marginal" effects.

- multinomial_probabilities:

  Logical; expand multinomial categories.

- show_coefficients:

  Logical; show edge coefficients.

## Value

A data frame suitable for ggplot2 rendering of a DAG.
