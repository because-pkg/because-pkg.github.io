# Sanitize Formula Term Name

Converts complex R formula terms into valid alphanumeric node/parameter
names for BUGS engines (JAGS/NIMBLE), parameter maps, d-separation
tests, and DAG layouts.

## Usage

``` r
sanitize_term_name(term)
```

## Arguments

- term:

  Character string (e.g., "A:B", "I(A^2)")

## Value

Sanitized string (e.g., "A_x_B", "A_pow2")
