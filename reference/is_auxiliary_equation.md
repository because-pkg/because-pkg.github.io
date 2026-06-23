# Create a Custom Covariance Structure for use with because

This helper function allows users to easily create custom covariance
structures by providing just a function that computes the precision (or
covariance) matrix. The function handles all the S3 method registration
automatically.

## Usage

``` r
is_auxiliary_equation(family, response, all_responses, ...)
```

## Arguments

- family:

  A family object (S3 class `because_family_*`)

- response:

  Character; the response variable name

- all_responses:

  Character vector of all response variable names in the model

- ...:

  Additional arguments

- name:

  Character string; the name of your structure (e.g., "spatial_knn").
  This will be used to create the S3 class and methods.

- precision_fn:

  A function that takes your structure data and returns a precision
  matrix (N x N). The first return value should be the precision matrix.
  Additional list elements can be returned and will be passed to JAGS as
  data.

- description:

  Optional description of the structure for documentation.

## Value

A constructor function that creates structure objects of your custom
class.

A list: `list(skip_likelihood = FALSE, skip_variance = FALSE)`

## Details

This function creates:

- A constructor function (returned) for creating structure objects

- S3 methods for `jags_structure_definition` (generates JAGS code)

- S3 methods for `prepare_structure_data` (prepares data for JAGS)

The precision function should return either:

- A single matrix (the precision matrix)

- A list with `Prec` (the precision matrix) and optionally other data
