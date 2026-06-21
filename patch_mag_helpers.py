import re

with open("R/mag_helpers.R", "r") as f:
    content = f.read()

# Replace the function signature
new_content = content.replace(
"""                                  categorical_vars = NULL,
                                  deterministic_terms = NULL,
                                  latent_sharers = list()) {""",
"""                                  categorical_vars = NULL,
                                  deterministic_terms = NULL,
                                  latent_sharers = list(),
                                  d_obj = NULL) {"""
)

# Insert the Ancestry Rule before the Root-node rule
ancestry_rule = """
        # Ancestry Rule: if one variable is an ancestor of the other in the DAG,
        # the descendant should be the Response (var1) and the ancestor should be the Predictor (var2).
        if (!is.null(d_obj)) {
            # dagitty::ancestors includes the node itself.
            if (var1 %in% names(d_obj) && var2 %in% names(d_obj)) {
                if (var1 %in% dagitty::ancestors(d_obj, var2)) {
                    # var1 is an ancestor of var2, so var2 is the descendant.
                    # Swap so var2 becomes Response (var1)
                    temp <- var1
                    var1 <- var2
                    var2 <- temp
                }
            }
        }
"""

new_content = new_content.replace(
    '        # Root-node rule: exogenous variables (no parents in the DAG) should',
    ancestry_rule + '\n        # Root-node rule: exogenous variables (no parents in the DAG) should'
)

with open("R/mag_helpers.R", "w") as f:
    f.write(new_content)

