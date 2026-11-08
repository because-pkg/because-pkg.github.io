# Sample from a compiled JAGS model

Sample from a compiled JAGS model

## Usage

``` r
sample_jags_model(model, monitor, n.iter, n.burnin, n.thin)
```

## Arguments

- model:

  Compiled jags.model object

- monitor:

  Character vector of variables to monitor

- n.iter:

  Number of iterations

- n.burnin:

  Number of burnin iterations

- n.thin:

  Thinning interval

## Value

mcmc.list of samples
