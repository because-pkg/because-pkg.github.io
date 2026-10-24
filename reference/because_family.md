# Create a Custom Family (Distribution) for use with because

This helper function allows users to easily create custom distribution
families by providing the JAGS likelihood code. The function handles S3
method registration.

## Usage

``` r
because_family(
  jags_likelihood = NULL,
  link = "identity",
  extra_priors = NULL,
  name = NULL,
  description = NULL
)
```

## Arguments

- jags_likelihood:

  A character string with the JAGS likelihood code. Use placeholders:
  `{response}` for variable name, `{mu}` for mean, `{i}` for loop index,
  `{tau}` for precision parameter.

- link:

  Link function name (default: "identity"). Options: "identity", "log",
  "logit".

- extra_priors:

  Character vector of additional JAGS prior statements. Use `{response}`
  placeholder for variable-specific naming.

- name:

  Optional character string; the name of your family (e.g.,
  "student_t"). If omitted, an internal name is generated automatically.

- description:

  Optional description of the family.

## Value

A constructor function for the family that can be passed to because(...,
family = ...).

## Details

**Engine compatibility:** Custom families defined via `because_family()`
inject BUGS-syntax likelihoods and priors into the model specification
and are supported for `engine = "jags"` and `engine = "nimble"`
(provided the BUGS distribution is available in NIMBLE). They are not
supported for `engine = "numpyro"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a Student-t family for robust regression
student_t <- because_family(
  name = "student_t",
  jags_likelihood = "{response}[{i}] ~ dt({mu}[{i}], {tau}, df_{response})",
  extra_priors = c("df_{response} ~ dunif(2, 100)")
)

# Use it
fit <- because(Y ~ X, data = data, family = student_t())
} # }
```
