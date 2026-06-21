import re

with open("R/because.R", "r") as f:
    content = f.read()

# Define the block to move
block_to_move = """
    # Auto-detect latent variables: variables in equations but not in data
    if (is.null(latent)) {
      vars_in_equations <- unique(unlist(lapply(equations, all.vars)))
      vars_in_data <- names(data)

      # Find variables that appear in equations but not in data
      potential_latents <- setdiff(vars_in_equations, vars_in_data)

      # Extension Hook: Remove specialized variables from potential latents
      potential_latents <- dsep_potential_latent_hook(
        family_obj,
        potential_latents
      )

      # Exclude polynomial internal variables (they're deterministic, not latent)
      if (!is.null(all_poly_terms)) {
        poly_internal_names <- sapply(all_poly_terms, function(x) {
          x$internal_name
        })
        potential_latents <- setdiff(potential_latents, poly_internal_names)
      }

      # Exclude categorical parent variables (they are replaced by dummies but are still in the model equations)
      cat_metadata <- attr(data, "categorical_vars")
      if (!is.null(cat_metadata)) {
        potential_latents <- setdiff(potential_latents, names(cat_metadata))
      }

      if (length(potential_latents) > 0) {
        # Auto-detect latent variables
        latent <- potential_latents

        if (!quiet) {
          message(
            "Auto-detected latent variable(s): ",
            paste(latent, collapse = ", "),
            "\\n(Variables in equations but not in data will be treated as latent.)\\n",
            "Generating m-separation tests for MAG..."
          )
        }
      }
    }

    if (!quiet) {
      if (!is.null(latent)) {
        message(
          "Generating m-separation tests (MAG with latent variables)..."
        )
      } else {
        message("Generating d-separation tests...")
      }
    }
"""

# The code string we'll replace might have slightly different indentation, so let's use a regex to capture it.
# Capture from "    # Auto-detect latent variables: variables in equations but not in data" up to before "# Expanded random terms for d-sep"

match = re.search(r'    # Auto-detect latent variables: variables in equations but not in data.*?    # Expanded random terms for d-sep', content, re.DOTALL)
if match:
    extracted_block = match.group(0).replace('    # Expanded random terms for d-sep', '')
    # Remove it from its original location
    new_content = content.replace(extracted_block, '\n')
    
    # Now inject it BEFORE "  # Handle d-sep logic"
    injection_point = '  # Handle d-sep logic'
    
    # We should also fix the message "Generating m-separation tests for MAG..." which is printed inside it,
    # and maybe only print "Generating ... " once when we actually do dsep?
    # Wait, the original code had those messages inside. Let's just keep them.
    
    new_content = new_content.replace(injection_point, extracted_block + '\n' + injection_point)
    
    with open("R/because.R", "w") as f:
        f.write(new_content)
    print("Success")
else:
    print("Could not find the block to extract!")
