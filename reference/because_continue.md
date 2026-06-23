# Continue MCMC Sampling for a Fitted Model

Extends an existing `because` model by running additional MCMC
iterations, then (optionally) combining the new samples with the
original ones. This is useful when traceplots or Rhat values suggest
that chains have not converged and more iterations are needed — without
discarding the burn-in and samples already obtained.

## Usage

``` r
because_continue(
  object,
  n.iter = 5000,
  combine = TRUE,
  parallel = NULL,
  n.cores = NULL,
  n.adapt = 200,
  quiet = FALSE
)
```

## Arguments

- object:

  A fitted `because` model (class `"because"`). Must have been run with
  `dsep = FALSE` (i.e. a full MCMC fit, not a d-separation-only run).

- n.iter:

  Integer. Number of *additional* post-warmup iterations to draw per
  chain, in the same units as the `n.iter` argument of
  [`because()`](https://because-pkg.github.io/because/reference/because.md).
  The number of stored samples added per chain is `n.iter %/% thin`,
  where `thin` is inherited from the original run.

- combine:

  Logical (default `TRUE`). If `TRUE` the new samples are appended to
  the existing ones so the returned object contains all iterations. If
  `FALSE` only the new samples are kept.

- parallel:

  Logical or `NULL` (default). `NULL` reuses the same setting as the
  original
  [`because()`](https://because-pkg.github.io/because/reference/because.md)
  call (stored in `object$parallel`). Set explicitly to `TRUE` or
  `FALSE` to override. Only relevant for JAGS rebuild paths; ignored for
  Nimble (see Details).

- n.cores:

  Integer or `NULL`. Number of cores for parallel execution. `NULL`
  reuses the value from the original run.

- n.adapt:

  Integer. Adaptation iterations used when the JAGS model must be
  rebuilt from scratch (ignored when the live model object is reused,
  and not applicable to Nimble). Default `200`.

- quiet:

  Logical. Suppress progress messages. Default `FALSE`.

## Value

The updated `because` object with extended (or replaced) samples and a
recalculated summary including Rhat. The live model object
(`result$model` for JAGS, `result$nimble_compiled` for Nimble) is
updated so that `because_continue()` can be called again without
recompilation.

## Details

**Engine support and practical guidance**

`because_continue()` is most useful with the **JAGS** engine, which is
the slowest of the three supported engines. The function provides a
*warm-start rebuild* path that works even across R session restarts,
making it a practical tool for extending long-running JAGS fits without
discarding the existing samples or redoing burn-in.

|  |  |  |
|----|----|----|
| **Engine** | **When it works** | **What happens** |
| `jags` | Always (see strategies below) | No or minimal burn-in |
| `nimble` | Same R session + `parallel = FALSE` only | `run(reset=FALSE)`, no recompilation |
| `numpyro` | Not supported | Re-run [`because()`](https://because-pkg.github.io/because/reference/because.md) with larger `n.iter` |

**JAGS — sampling strategies (tried in order)**

1.  **Live model (sequential runs, same session)** — if `object$model`
    is still a valid `jags` object,
    [`rjags::coda.samples()`](https://rdrr.io/pkg/rjags/man/coda.samples.html)
    is called directly. No recompilation, no burn-in, exact continuation
    of all chains.

2.  **Parallel rebuild** — if `parallel = TRUE`, each chain is
    recompiled in its own worker process and warm-started from the last
    sample of that chain (`alpha_*` and `beta_*` parameters only, since
    derived nodes such as `sigma_*_res` or `lambda_*` cannot be
    initialised directly in JAGS).

3.  **Sequential rebuild** — a single `jags.model()` call with
    warm-start inits. Automatically falls back to progressively safer
    init sets (alpha/beta only, then cold start) when deterministic
    nodes are encountered.

The rebuild paths work *across R sessions* because the JAGS model code
is stored in `object$model_code` and recompilation is fast (seconds, not
minutes).

**Nimble — live compiled object path**

For Nimble models compiled in the *current* R session with
`parallel = FALSE`, the C++ compiled MCMC object
(`object$nimble_compiled`) is reused directly via
`compiled_mcmc$run(reset = FALSE)`. Each chain is reinitialised to the
last sample of that chain in the stored `mcmc.list`, so the continuation
starts from the exact posterior position — no recompilation (which can
take several minutes for large models).

Nimble continuation is **not available** when:

- The original run used `parallel = TRUE` — compiled objects live in
  worker processes and are destroyed when sampling finishes.

- R has been restarted — C++ compiled objects are session-specific and
  cannot be saved to disk.

In those cases, re-run
[`because()`](https://because-pkg.github.io/because/reference/because.md)
with a larger `n.iter`.

**NumPyro**

NumPyro uses HMC/NUTS, which converges so efficiently that extending an
existing run is rarely necessary. Re-run
[`because()`](https://because-pkg.github.io/because/reference/because.md)
with a larger `n.iter` if more samples are needed.

## Examples

``` r
if (FALSE) { # \dontrun{
# --- JAGS: sequential run, same session ---
# Live model object is reused — instant continuation, no recompilation.
fit <- because(equations = list(y ~ x), data = dat,
               n.iter = 2000, n.burnin = 500, n.chains = 3)
summary(fit)  # Rhat looks borderline
fit <- because_continue(fit, n.iter = 5000)
summary(fit)  # now 6500 effective iterations

# --- JAGS: after restarting R ---
# Model is rebuilt from object$model_code with warm-start inits.
# parallel = NULL reuses the original setting automatically.
load("fit.RData")
fit <- because_continue(fit, n.iter = 5000)

# --- Nimble: sequential run, same session ---
fit_nim <- because(equations = list(y ~ x), data = dat,
                   engine = "nimble", n.iter = 2000, n.burnin = 500,
                   parallel = FALSE)
fit_nim <- because_continue(fit_nim, n.iter = 5000)
} # }
```
