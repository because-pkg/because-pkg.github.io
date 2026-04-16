# Generate a JAGS model string for Bayesian SEM (Because)

This function builds the model code to be passed to JAGS based on a set
of structural equations. It supports custom covariance structures
(spatial, phylogenetic, etc.). Missing values are handled both in the
response and predictor variables treating all of them as stochastic
nodes.

## Usage

``` r
because_model(
  equations,
  is_multi_structure = FALSE,
  latent_method = "correlations",
  structures = list(),
  random_structure_names = NULL,
  random_terms = list(),
  vars_with_na = NULL,
  induced_correlations = NULL,
  variability = NULL,
  family = NULL,
  standardize_latent = TRUE,
  poly_terms = NULL,
  latent = NULL,
  categorical_vars = NULL,
  fix_latent = "loading",
  fix_residual_variance = NULL,
  priors = NULL,
  hierarchical_info = NULL,
  engine = "jags",
  quiet = FALSE
)
```

## Arguments

- equations:

  A list of model formulas.

- is_multi_structure:

  Logical (Internal). If TRUE, handles 3D Precision arrays.

- latent_method:

  Method for handling latent variables ("correlations" or "explicit").

- structures:

  A named list of structural objects (e.g. matrices, trees) to include
  as correlations.

- random_structure_names:

  Optional character vector of structural names applied to all
  variables.

- random_terms:

  Optional list of random effects (response, group).

- vars_with_na:

  Optional character vector of response variable names that have missing
  data. These variables will use element-wise likelihoods instead of
  multivariate normal.

- induced_correlations:

  Optional list of variable pairs with induced correlations from latent
  variables. Each element should be a character vector of length 2
  specifying the pair of variables that share a latent common cause.

- variability:

  Optional character vector or named character vector of variable names
  that have measurement error or within-species variability. If named,
  the names should be the variable names and the values should be the
  type of variability: "se" (for mean and standard error) or "reps" (for
  repeated measures). If unnamed, it defaults to "se" for all specified
  variables.

  - "se": Expects `Var_mean` and `Var_se` in the data. The model fixes
    observation error: `Var_mean ~ dnorm(Var, 1/Var_se^2)`.

  - "reps": Expects `Var_obs` (matrix) and `N_reps_Var` (vector) in the
    data. The model estimates observation error:
    `Var_obs[i,j] ~ dnorm(Var[i], Var_tau)`.

- family:

  Optional named character vector specifying the family/distribution for
  response variables. Default is "gaussian" for all variables. Supported
  values: "gaussian", "binomial", "multinomial". For "binomial"
  variables, the model uses a logit link and a Bernoulli likelihood,
  with phylogenetic correlation modeled on the latent scale.

- standardize_latent:

  Logical (default TRUE). If TRUE, standardizes latent variables to unit
  variance.

- poly_terms:

  (Internal) List of polynomial terms for model generation.

- latent:

  Optional character vector of latent variable names.

- categorical_vars:

  Optional character vector of categorical variable names.

- fix_residual_variance:

  Optional numeric value or named vector to fix residual variance.

- priors:

  Optional named list of custom priors.

- hierarchical_info:

  Optional list containing data hierarchy (levels, link_vars).

- engine:

  Bayesian engine to use ("jags" or "nimble").

## Value

A list with two elements:

- `model`: A character string containing the JAGS model code.

- `parameter_map`: A data frame mapping response variables to their
  predictors and parameter names.

## Details

The generated model includes:

- Linear predictors and multivariate normal likelihoods for each
  response variable.

- Priors for intercepts (`alpha`), slopes (`beta`), and residual
  precisions (`tau`).

- Custom covariance modeled via provided structural objects (e.g. VCV
  matrices).

- (Optional) Observation models for variables with measurement error:

  - Type "se": `Var_mean ~ dnorm(Var, 1/Var_se^2)`

  - Type "reps": `Var_obs[i,j] ~ dnorm(Var[i], Var_tau)`

- (Optional) Generalized linear mixed models for non-Gaussian responses
  (e.g., binomial).

- (Optional) Element-wise likelihoods for response variables with
  missing data.

## Examples

``` r
eqs <- list(BR ~ BM, S ~ BR, G ~ BR, L ~ BR)
cat(because_model(eqs, is_multi_structure = TRUE)$model)
#> model {
#>   # Common structures and priors
#>   # Structural equations
#> 
#>   for (i in 1:N) {
#>     mu_BR[i] <- alpha_BR + beta_BR_BM*BM[i]
#>   }
#>   for (i in 1:N) {
#>     mu_S[i] <- alpha_S + beta_S_BR*BR[i]
#>   }
#>   for (i in 1:N) {
#>     mu_G[i] <- alpha_G + beta_G_BR*BR[i]
#>   }
#>   for (i in 1:N) {
#>     mu_L[i] <- alpha_L + beta_L_BR*BR[i]
#>   }
#>   # Multivariate normal likelihoods
#>   for (i in 1:N) {
#>     BR[i] ~ dnorm(mu_BR[i], tau_res_BR)
#>     log_lik_BR[i] <- logdensity.norm(BR[i], mu_BR[i], tau_res_BR)
#>   }
#>   for (i in 1:N) {
#>     S[i] ~ dnorm(mu_S[i], tau_res_S)
#>     log_lik_S[i] <- logdensity.norm(S[i], mu_S[i], tau_res_S)
#>   }
#>   for (i in 1:N) {
#>     G[i] ~ dnorm(mu_G[i], tau_res_G)
#>     log_lik_G[i] <- logdensity.norm(G[i], mu_G[i], tau_res_G)
#>   }
#>   for (i in 1:N) {
#>     L[i] ~ dnorm(mu_L[i], tau_res_L)
#>     log_lik_L[i] <- logdensity.norm(L[i], mu_L[i], tau_res_L)
#>   }
#>   # Priors for structural parameters
#>   alpha_BR ~ dnorm(0, 0.01)
#>   sigma_BR_res ~ dunif(0, 100)
#>   tau_res_BR <- 1 / (sigma_BR_res * sigma_BR_res)
#>   alpha_S ~ dnorm(0, 0.01)
#>   sigma_S_res ~ dunif(0, 100)
#>   tau_res_S <- 1 / (sigma_S_res * sigma_S_res)
#>   alpha_G ~ dnorm(0, 0.01)
#>   sigma_G_res ~ dunif(0, 100)
#>   tau_res_G <- 1 / (sigma_G_res * sigma_G_res)
#>   alpha_L ~ dnorm(0, 0.01)
#>   sigma_L_res ~ dunif(0, 100)
#>   tau_res_L <- 1 / (sigma_L_res * sigma_L_res)
#>   beta_BR_BM ~ dnorm(0, 1.0E-6)
#>   beta_S_BR ~ dnorm(0, 1.0E-6)
#>   beta_G_BR ~ dnorm(0, 1.0E-6)
#>   beta_L_BR ~ dnorm(0, 1.0E-6)
#>   # Phylogenetic uncertainty weighting
#>   for (k in 1:Ntree) {
#>     p_tree[k] <- 1/Ntree
#>   }
#>   K ~ dcat(p_tree[1:Ntree])
#> }
```
