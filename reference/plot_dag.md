# Plot DAG from Equations or Fitted Model

Visualizes the Directed Acyclic Graph (DAG) implied by a set of
equations or a fitted `because` model. If a fitted model is provided:

- Path coefficients are displayed on the edges.

- Edge colors are determined by the `edge_color_scheme` argument (e.g.,
  Red/Blue for directional, Black/Grey for binary).

- Structural edges (without fitted data) are colored black.

- Edge thickness scales with the absolute effect size.

## Usage

``` r
plot_dag(
  x,
  layout = "kk",
  latent = NULL,
  node_size = 14,
  node_color = "black",
  node_fill = "white",
  node_stroke = 1.5,
  text_size = 4,
  edge_width_range = c(0.5, 2),
  edge_color_scheme = c("directional", "binary", "monochrome"),
  show_coefficients = TRUE,
  coords = NULL,
  family = NULL
)
```

## Arguments

- x:

  A list of formulas (equations), a `because` model object, or a list of
  these.

- layout:

  The layout algorithm to use for positioning nodes (default `"kk"`,
  Kamada-Kawai). Any layout name accepted by
  [`create_layout`](https://ggraph.data-imaginist.com/reference/ggraph.html)
  can be used, including:

  - `"kk"` — Kamada-Kawai spring layout (default; good general-purpose
    layout).

  - `"sugiyama"` — hierarchical/layered layout; useful for simple
    chains.

  - `"fr"` — Fruchterman-Reingold spring layout.

  - `"nicely"` — automatic choice based on graph properties.

  - `"circle"` — nodes arranged in a circle.

  Override node positions entirely with the `coords` argument.

- latent:

  Character vector of latent variable names. Overrides the model's
  latent variables if provided.

- node_size:

  Size of the nodes (default 14).

- node_color:

  Color of the node border (default "black").

- node_fill:

  Color of the node interior (default "white").

- node_stroke:

  Thickness of the node border (default 1.5).

- text_size:

  Size of the labels (default 4).

- edge_width_range:

  Vector of length 2 defining the range of arrow widths (min, max) based
  on effect size.

- edge_color_scheme:

  Character; one of "directional" (default), "binary", or "monochrome".
  "directional" colors edges red/blue/grey based on effect direction and
  whether the 95\\ "binary" colors edges black/grey based on whether the
  95\\ "monochrome" colors all edges black.

- show_coefficients:

  Logical; whether to print coefficient values on edges (only for fitted
  models).

- coords:

  Optional named list of coordinates for the nodes, e.g.
  `list(A = c(1, 1), B = c(2, 2))`. If provided, these will override the
  `layout` algorithm.

- family:

  Optional named character vector of families for response variables.

## Value

A `ggplot` object that can be further customized with standard ggplot2
functions (e.g., `+ theme_...()`, `+ ggtitle(...)`).

## Details

Interaction terms (e.g. `BM:M`) and
[`I()`](https://rdrr.io/r/base/AsIs.html) transformations are rendered
as **explicit intermediate nodes** (grey diamonds), following the
Interaction DAG (IDAG) convention of Attia, Holliday & Oldmeadow (2022).
This makes the deterministic nature of these terms visually clear and is
consistent with how `because_dsep` treats them for d-separation.

## References

Attia, J., Holliday, E., & Oldmeadow, C. (2022). A proposal for
capturing interaction and effect modification using DAGs. *International
Journal of Epidemiology*, 51(4), 1047–1053.
[doi:10.1093/ije/dyac105](https://doi.org/10.1093/ije/dyac105)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic plotting
eq <- list(y ~ x + z, x ~ z)
plot_dag(eq)

# Custom Layout
my_coords <- list(
  y = c(1, 1),
  x = c(0, 0),
  z = c(2, 0)
)
plot_dag(eq, coords = my_coords)
} # }
```
