            # Independent Model (No random effects)
            model_lines <- c(
              model_lines,
              paste0("  for (i in 1:", loop_bound, ") {"),
              paste0(
                "    ",
                response_var,
                "[i] ~ dnorm(",
                mu,
                "[i], ",
                tau_e,
                ")"
              )
            )
            if (engine == "jags") {
              model_lines <- c(
                model_lines,
                paste0(
                  "    log_lik_",
                  response,
                  suffix,
                  "[i] <- logdensity.norm(",
                  response_var,
                  "[i], ",
                  mu,
                  "[i], ",
                  tau_e,
                  ")"
                )
              )
            }
            model_lines <- c(model_lines, "  }")
          } else {
            # Optimized Random Effects Formulation (Additive)
            model_lines <- c(
              model_lines,
              paste0("  for (i in 1:", loop_bound, ") {")
            )
            additive_terms <- ""
