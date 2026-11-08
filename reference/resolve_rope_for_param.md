# Resolve ROPE interval for a specific parameter

Resolves in priority:

1.  Exact parameter match (e.g. rope\["beta_abundance_prey"\])

2.  Response variable match (e.g. rope\["abundance"\])

3.  Distribution family match (e.g. rope\["poisson"\])

4.  List default (rope\["default"\])

5.  Principled family default (get_family_default_rope)

## Usage

``` r
resolve_rope_for_param(pname, resp = NULL, fam = NULL, rope_spec = NULL)
```

## Arguments

- pname:

  Parameter name

- resp:

  Response variable name

- fam:

  Family name

- rope_spec:

  User-supplied ROPE specification (NULL, TRUE, vector, or list)

## Value

Numeric vector of length 2
