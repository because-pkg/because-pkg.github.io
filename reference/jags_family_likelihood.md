# Generate BUGS likelihood code for a distribution family

This generic allows S3 dispatch to generate the appropriate BUGS
likelihood code (for JAGS and NIMBLE) for any distribution family.

## Usage

``` r
jags_family_likelihood(
  family,
  response,
  predictors = NULL,
  suffix = "",
  has_structure = FALSE,
  link = "identity",
  ...
)
```

## Arguments

- family:

  A family object (created by
  [`get_family_object`](https://because-pkg.github.io/because/reference/get_family_object.md)
  or custom constructor)

- response:

  Character string; name of the response variable

- predictors:

  Character vector; names of predictor variables (may be NULL)

- suffix:

  Character string; suffix for variable names (e.g., "1" for multiple
  responses)

- has_structure:

  Logical; whether the model includes a structure (e.g., phylogenetic)

- link:

  Character string; link function ("identity", "log", "logit")

- ...:

  Additional arguments passed to methods

## Value

A list with:

- likelihood_code:

  Character vector of BUGS likelihood statements

- prior_code:

  Character vector of BUGS prior statements (or NULL)

- data_requirements:

  Character vector of required data elements (or NULL)
