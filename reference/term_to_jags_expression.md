# Convert R Term to JAGS Expression

Transforms R syntax into JAGS-compatible deterministic code

## Usage

``` r
term_to_jags_expression(term)
```

## Arguments

- term:

  Original R term (e.g., "A:B")

## Value

JAGS expression (e.g., "A\\i\] \* B\\i\]")
