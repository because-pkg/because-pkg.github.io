# Compile and run a NIMBLE model with standard because settings

Compile and run a NIMBLE model with standard because settings

## Usage

``` r
run_nimble_model(
  nimble_code,
  constants,
  data,
  inits,
  n.chains,
  n.iter,
  n.burnin,
  n.thin,
  monitor,
  quiet,
  family = NULL,
  nimble_samplers = NULL,
  nimble_waic = TRUE
)
```

## Arguments

- nimble_code:

  NIMBLE model code (from nimbleCode)

- constants:

  Data constants for NIMBLE

- data:

  NIMBLE data list

- inits:

  Initial values list

- n.chains:

  Number of MCMC chains

- n.iter:

  Number of iterations

- n.burnin:

  Number of burnin iterations

- n.thin:

  Thinning interval

- monitor:

  Character vector of parameters to monitor

- quiet:

  Logical, suppress output

- family:

  Named vector of family types for hardening samplers

- nimble_samplers:

  Optional custom sampler list

- nimble_waic:

  Logical, calculate WAIC

## Value

A list containing samples and the compiled model
