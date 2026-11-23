# Convert NumPyro output into coda mcmc.list

Filters sampled NumPyro arrays into a
[`coda::mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html)
according to `monitor` ("interpretable", "all", or explicit parameter
names), retaining all structural vectors regardless of length and
informing the user when large random-effect or latent vectors are
omitted for memory safety.

## Usage

``` r
format_numpyro_samples(py_result, monitor = "interpretable", quiet = FALSE)
```

## Arguments

- py_result:

  Result list from fit_numpyro_model

- monitor:

  Monitoring mode ("interpretable", "all", or character vector of
  parameter names)

- quiet:

  Logical; if FALSE, emits an informational message when large vector
  parameters (e.g. high-dimensional `u_*` or imputed node vectors) are
  omitted under `monitor = "interpretable"`.

## Value

coda::mcmc.list
