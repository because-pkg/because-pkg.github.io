# Execute full JAGS pipeline (sequential or parallel) for because()

Execute full JAGS pipeline (sequential or parallel) for because()

## Usage

``` r
run_jags_pipeline(
  model_file,
  model_string,
  data,
  extension_inits,
  monitor,
  n.chains,
  n.iter,
  n.burnin,
  n.thin,
  n.adapt,
  quiet,
  verbose,
  parallel,
  n.cores,
  cl = NULL,
  DIC = FALSE,
  WAIC = FALSE
)
```
