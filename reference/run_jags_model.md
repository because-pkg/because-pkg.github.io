# Run a JAGS model with standard because settings

Run a JAGS model with standard because settings

## Usage

``` r
run_jags_model(
  model_file,
  data,
  inits_list,
  n.chains,
  n.adapt,
  quiet,
  model_string
)
```

## Arguments

- model_file:

  Path to the JAGS model file or text connection

- data:

  List of data for the model

- inits_list:

  List of initial values

- n.chains:

  Number of MCMC chains

- n.adapt:

  Number of adaptation steps

- quiet:

  Logical, suppress output

- model_string:

  The model string (for error reporting)

## Value

A compiled jags.model object
