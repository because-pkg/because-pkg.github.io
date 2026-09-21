#' Run d-separation tests or equation adjustments for a because model
#'
#' Orchestrates d-separation testing for NumPyro or JAGS/NIMBLE backends,
#' or performs equation adjustments (MAG marginalization and categorical expansion)
#' when dsep = FALSE.
#'
#' @keywords internal
run_because_dsep <- function(
  dsep, engine, equations, data, family, structures, structure_obj,
  hierarchical_info, is_hierarchical, random_terms, levels, multiscale,
  link_vars, latent, latent_method, id_col, variability, all_poly_terms,
  fixed_equations_temp, induced_cors_in = NULL, dsep_max_obs, aggregate_crossscale,
  parallel, n.cores, n.chains, n.iter, n.burnin, n.thin, n.adapt,
  ic_recompile, fix_residual_variance, quiet, priors, monitor_mode,
  expand_ordered, nimble_samplers, adapt_delta, max_treedepth,
  prior_scale_fixed, verbose,
  # Additional variables from parent scope
  family_obj = NULL, hierarchy = NULL, original_call = NULL,
  original_data = NULL, random = NULL, response_vars_with_na = NULL,
  reuse_models = FALSE, structure = NULL, because_py = NULL,
  WAIC = FALSE, DIC = FALSE
) {
  if (dsep) {
    if (engine == "numpyro") {
      return(run_dsep_numpyro(
        equations = equations,
        structures = structures,
        data = data,
        latent = latent,
        random_terms = random_terms,
        hierarchical_info = hierarchical_info,
        all_poly_terms = all_poly_terms,
        family = family,
        quiet = quiet,
        original_call = original_call,
        original_data = original_data,
        id_col = id_col,
        variability = variability,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        expand_ordered = expand_ordered,
        reuse_models = reuse_models,
        aggregate_crossscale = aggregate_crossscale,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.chains = n.chains,
        n.thin = n.thin,
        parallel = parallel,
        n.cores = n.cores,
        dsep_max_obs = dsep_max_obs,
        prior_scale_fixed = prior_scale_fixed,
        engine = engine,
        structure = structure,
        because_py = because_py
      ))
    } else {
      return(run_dsep_jags_nimble(
        engine = engine,
        equations = equations,
        data = data,
        family = family,
        structures = structures,
        structure_obj = structure_obj,
        hierarchical_info = hierarchical_info,
        is_hierarchical = is_hierarchical,
        random_terms = random_terms,
        levels = levels,
        multiscale = multiscale,
        link_vars = link_vars,
        latent = latent,
        latent_method = latent_method,
        id_col = id_col,
        variability = variability,
        all_poly_terms = all_poly_terms,
        fixed_equations_temp = fixed_equations_temp,
        induced_cors_in = induced_cors_in,
        dsep_max_obs = dsep_max_obs,
        aggregate_crossscale = aggregate_crossscale,
        parallel = parallel,
        n.cores = n.cores,
        n.chains = n.chains,
        n.iter = n.iter,
        n.burnin = n.burnin,
        n.thin = n.thin,
        n.adapt = n.adapt,
        ic_recompile = ic_recompile,
        fix_residual_variance = fix_residual_variance,
        quiet = quiet,
        priors = priors,
        monitor_mode = monitor_mode,
        expand_ordered = expand_ordered,
        nimble_samplers = nimble_samplers,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        prior_scale_fixed = prior_scale_fixed,
        verbose = verbose,
        family_obj = family_obj,
        hierarchy = hierarchy,
        original_call = original_call,
        original_data = original_data,
        random = random,
        response_vars_with_na = response_vars_with_na,
        reuse_models = reuse_models,
        structure = structure,
        WAIC = WAIC,
        DIC = DIC
      ))
    }
  }

  # When dsep = FALSE: Handle latent variable method (MAG vs explicit)
  mag_res <- process_latent_mag_equations(
    latent = latent,
    latent_method = latent_method,
    equations = equations,
    random_terms = random_terms,
    hierarchical_info = hierarchical_info,
    family = family,
    quiet = quiet,
    induced_cors = induced_cors_in,
    dsep = FALSE
  )
  equations    <- mag_res$equations
  induced_cors <- mag_res$induced_cors

  # Auto-expand categorical variables in equations
  cat_res <- expand_categorical_equations(
    equations = equations,
    data = data,
    response_vars_with_na = response_vars_with_na,
    quiet = quiet
  )
  equations <- cat_res$equations
  data      <- cat_res$data

  return(list(
    dsep_tests   = NULL,
    induced_cors = induced_cors,
    equations    = equations,
    data         = data
  ))
}
