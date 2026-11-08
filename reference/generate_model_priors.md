# Generate priors and variance partitioning for JAGS/BUGS model

Constructs priors for intercepts (alpha), slopes (beta), residual
precisions (tau), structural components (tau_u, sigma), structural
variance partitioning (Variance Partition Coefficient / VPC, e.g.
Pagel's lambda in phylogenetic models or spatial VPC), and
family-specific auxiliary parameters.

## Usage

``` r
generate_model_priors(ctx)
```

## Arguments

- ctx:

  Builder context environment
