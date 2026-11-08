# Expand Categorical Variables in Structural Equations

Automatically expands categorical variables in model formulas into dummy
indicator variables, supporting interaction terms and deterministic
identity definitions for imputed factors.

## Usage

``` r
expand_categorical_equations(
  equations,
  data,
  response_vars_with_na = NULL,
  quiet = FALSE
)
```
