# Extract d-separation statements from a structural equation model

This function takes a set of structural equations defining a causal
model and returns the conditional independence statements (d-separation
or m-separation tests) implied by the model structure. If latent
variables are specified, the function uses the MAG (Mixed Acyclic Graph)
approach by Shipley and Douma (2021) to account for unmeasured latent
variables.

## Usage

``` r
because_dsep(
  equations,
  latent = NULL,
  random_terms = list(),
  hierarchical_info = NULL,
  poly_terms = NULL,
  categorical_vars = NULL,
  family = NULL,
  quiet = FALSE
)
```

## Arguments

- equations:

  A list of model formulas (one per structural equation), e.g.,
  `list(Y ~ X1 + X2, Z ~ Y)`.

- latent:

  Optional character vector of latent (unmeasured) variable names. If
  provided, the function converts the DAG to a MAG and returns
  m-separation tests.

- random_terms:

  Optional list of random effects (group, type) parsed from equations.

- hierarchical_info:

  Internal argument used to pass data hierarchy information (levels,
  grouping variables) for future implementation of multilevel
  d-separation tests (following Shipley 2009). Currently unused by the
  d-separation logic.

- quiet:

  Logical; if FALSE (default), print the basis set and MAG structure. If
  TRUE, suppress informational output.

## Value

If `latent` is NULL, returns a list of formulas representing conditional
independence tests. If `latent` is specified, returns a list with:

- `tests`: List of m-separation test formulas

- `correlations`: List of variable pairs with induced correlations

## Details

The function implements the basis set approach to d-separation testing
(Shipley 2000, 2009, 2016). For standard DAGs without latent variables,
it identifies pairs of non-adjacent variables and creates conditional
independence tests.

When latent variables are specified, the function uses the DAG-to-MAG
conversion (Shipley & Douma 2021) to identify m-separation statements
and induced correlations among observed variables that arise from shared
latent common causes.

Deterministic nodes (interaction terms such as `A:B`, and arithmetic
transformations such as `I(A^2)`) are kept as **explicit intermediate
nodes** in the DAG, following the D-separation extension of Geiger,
Verma & Pearl (1990). This ensures that the basis set includes
independence tests that condition on the deterministic term itself (e.g.
\\TL \perp BM \mid \\BM{:}M\\\\), which would be silently dropped if the
interaction were collapsed to its component variables.

## References

Geiger, D., Verma, T., & Pearl, J. (1990). Identifying independence in
Bayesian Networks. *Networks*, 20(5), 507–534.

Shipley, B. (2000). A new inferential test for path models based on
directed acyclic graphs. Structural Equation Modeling, 7(2), 206-218.

Shipley, B. (2009). Confirmatory path analysis in a generalized
multilevel context. Ecology, 90(2), 363-368.

Shipley, B. (2016). Cause and Correlation in Biology (2nd ed.).
Cambridge University Press.

Shipley, B., & Douma, J. C. (2021). Testing Piecewise Structural
Equations Models in the Presence of Latent Variables and Including
Correlated Errors. Structural Equation Modeling: A Multidisciplinary
Journal, 28(4), 582-589. https://doi.org/10.1080/10705511.2020.1871355

## Examples

``` r
# Standard DAG
equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)
ind_tests <- because_dsep(equations)
#> Basis Set for DAG: 
#> I(X,Y|Z) means X is d-separated from Y given the set Z in the DAG 
#> I( RS , BM |  ) 
#> I( DD , RS | NL ) 
#> I( LS , RS | BM ) 
#> I( DD , BM | NL ) 
#> I( NL , LS | BM, RS ) 
#> I( DD , LS | BM, NL ) 

# With latent variable
equations_latent <- list(X ~ Quality, Y ~ Quality)
result <- because_dsep(equations_latent, latent = "Quality")
#> Basis Set for MAG: 
#> I(X,Y|Z) means X is m-separated from Y given the set Z in the MAG 
#> No elements in the basis set 
# result$tests: m-separation tests
# result$correlations: induced correlation between X and Y
```
