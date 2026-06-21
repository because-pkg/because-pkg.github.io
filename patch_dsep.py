with open("R/because_dsep.R", "r") as f:
    content = f.read()

# For dsep_standard:
content = content.replace(
"""  tests <- mag_basis_to_formulas(
    basis,
    categorical_vars = categorical_vars,
    family = family,
    deterministic_terms = all_det_terms,
    root_vars = root_vars,
    hierarchical_info = hierarchical_info
  )""",
"""  tests <- mag_basis_to_formulas(
    basis,
    categorical_vars = categorical_vars,
    family = family,
    deterministic_terms = all_det_terms,
    root_vars = root_vars,
    hierarchical_info = hierarchical_info,
    d_obj = d_obj
  )""")

# For because_dsep (with latents):
content = content.replace(
"""  tests <- mag_basis_to_formulas(
    basis,
    latent_children = latent_children,
    categorical_vars = categorical_vars,
    family = family,
    deterministic_terms = all_det_terms,
    root_vars = root_vars,
    hierarchical_info = hierarchical_info,
    latent_sharers = latent_sharers
  )""",
"""  tests <- mag_basis_to_formulas(
    basis,
    latent_children = latent_children,
    categorical_vars = categorical_vars,
    family = family,
    deterministic_terms = all_det_terms,
    root_vars = root_vars,
    hierarchical_info = hierarchical_info,
    latent_sharers = latent_sharers,
    d_obj = d_obj
  )""")

with open("R/because_dsep.R", "w") as f:
    f.write(content)
