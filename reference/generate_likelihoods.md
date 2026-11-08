# Generate observation likelihoods for JAGS/BUGS model

Constructs the observation likelihood statements for all response
variables, including structured errors (phylogenetic, spatial,
multiPhylo), unstructured random intercepts, and MAG correlated errors.

## Usage

``` r
generate_likelihoods(ctx)
```

## Arguments

- ctx:

  Builder context environment

## Details

Dispatches to per-family helpers defined in:

- likelihood_continuous.R (gaussian)

- likelihood_categorical.R (multinomial, ordinal)

- likelihood_discrete.R (binomial, poisson, negbinomial/zinb,
  extensions)
