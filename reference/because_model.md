# Generate BUGS Model Code and Parameter Map for Bayesian SEM

Builds the BUGS-dialect model specification (used directly by JAGS and
parsed into an R expression by NIMBLE) along with the canonical
`parameter_map` data frame used across all engines (JAGS, NIMBLE, and
NumPyro). Supports custom covariance structures (spatial, phylogenetic,
etc.) and handles missing values in both response and predictor
variables as stochastic nodes.

## Usage

``` r
because_model(
  equations,
  is_multi_structure = FALSE,
  latent_method = "correlations",
  structures = list(),
  random_structure_names = NULL,
  random_terms = list(),
  vars_with_na = NULL,
  induced_correlations = NULL,
  variability = NULL,
  family = NULL,
  standardize_latent = TRUE,
  poly_terms = NULL,
  latent = NULL,
  categorical_vars = NULL,
  fix_latent = "loading",
  fix_residual_variance = NULL,
  priors = NULL,
  hierarchical_info = NULL,
  engine = "jags",
  quiet = FALSE
)
```

## Arguments

- equations:

  A list of model formulas.

- is_multi_structure:

  Logical; if `TRUE`, handles 3D Precision arrays (e.g. multi-object
  sampling).

- latent_method:

  Method for handling latent variables ("correlations" or "explicit").

- structures:

  A named list of structural objects (e.g. matrices, trees) to include
  as correlations.

- random_structure_names:

  Optional character vector of structural names applied to all
  variables.

- random_terms:

  Optional list of random effects (response, group).

- vars_with_na:

  Optional character vector of response variable names that have missing
  data.

- induced_correlations:

  Optional list of variable pairs with induced correlations from latent
  variables.

- variability:

  Optional character vector or named character vector of variable names
  that have measurement error or within-unit variability.

- family:

  Optional named character vector specifying the family/distribution for
  response variables.

- standardize_latent:

  Logical (default TRUE). If TRUE, standardizes latent variables to unit
  variance.

- poly_terms:

  (Internal) List of polynomial terms for model generation.

- latent:

  Optional character vector of latent variable names.

- categorical_vars:

  Optional character vector of categorical variable names.

- fix_latent:

  Variable to fix for latent identification.

- fix_residual_variance:

  Optional numeric value or named vector to fix residual variance.

- priors:

  Optional named list of custom priors.

- hierarchical_info:

  Optional list containing multiscale data hierarchy (levels,
  link_vars).

- engine:

  Bayesian engine to use (`"jags"`, `"nimble"`, or `"numpyro"`).

- quiet:

  Logical; suppress status messages.

## Value

A list with two elements:

- `model`: A character string containing the BUGS/JAGS model code.

- `parameter_map`: A data frame mapping response variables to their
  predictors and parameter names.
