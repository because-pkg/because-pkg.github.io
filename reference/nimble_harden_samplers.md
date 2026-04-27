# Harden NIMBLE Sampler Configuration

Applies robust sampler assignments to a NIMBLE MCMC configuration
object. Exported as a standalone function so that parallel worker nodes
— which load the package fresh — always use the current installed
version of this logic, regardless of which version of
[`because()`](https://because-pkg.github.io/because/reference/because.md)
was originally called.

## Usage

``` r
nimble_harden_samplers(
  mcmc_conf,
  family = NULL,
  nimble_samplers = NULL,
  quiet = TRUE
)
```

## Arguments

- mcmc_conf:

  A NIMBLE MCMC configuration object (from `configureMCMC()`).

- family:

  Named character vector of response families (same as
  [`because()`](https://because-pkg.github.io/because/reference/because.md)).

- nimble_samplers:

  Optional named list of user-specified samplers.

- quiet:

  Logical. If TRUE, suppress status messages.

## Value

The modified `mcmc_conf` object (invisibly).
