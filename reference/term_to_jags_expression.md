# Convert R Formula Term to BUGS Expression

Transforms R formula syntax into indexed BUGS deterministic expressions
for JAGS and NIMBLE.

## Usage

``` r
term_to_jags_expression(term)
```

## Arguments

- term:

  Original R term (e.g., "A:B")

## Value

BUGS expression (e.g., "A\\i\] \* B\\i\]")
