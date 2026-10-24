#' Prepare MCMC Monitors, Inits, and Clean Data
#'
#' Resolves monitored parameters according to mode ("interpretable", "all", or custom),
#' adds WAIC log_lik parameters if requested, stabilizes NIMBLE precision matrices,
#' extracts family-specific initial values, and prunes unused data fields.
#'
#' @noRd
prepare_monitors_and_inits <- function(
  model_string,
  monitor,
  response_vars,
  mag_exogenous_vars = character(0),
  family_obj = NULL,
  equations = list(),
  engine = "jags",
  data = list(),
  variability = NULL,
  variability_list = list(),
  WAIC = FALSE,
  quiet = FALSE
) {
  # Handle monitor mode
  monitor_mode <- NULL
  custom_monitors <- character(0)

  if (!is.null(monitor)) {
    if (is.character(monitor)) {
      if ("interpretable" %in% monitor) {
        monitor_mode <- "interpretable"
        custom_monitors <- setdiff(monitor, "interpretable")
      } else if ("all" %in% monitor) {
        monitor_mode <- "all"
        custom_monitors <- setdiff(monitor, "all")
      } else if (identical(monitor, "")) {
        monitor_mode <- "interpretable"
      } else {
        # Entirely custom vector
        custom_monitors <- monitor
      }
    }
  } else {
    monitor_mode <- "interpretable"
  }

  lines <- unlist(strsplit(model_string, "\n"))

  extract_names <- function(pattern) {
    out <- grep(pattern, lines, value = TRUE)
    out <- grep("(<-|~)", out, value = TRUE)
    matches <- regmatches(
      out,
      regexec(
        "(?:logit|log|cloglog|probit)?\\(?\\s*([a-zA-Z0-9_.]+)(?:\\[.*\\])?\\)?\\s*(?:<-|~)",
        out
      )
    )

    res <- sapply(matches, function(m) {
      if (length(m) >= 2) m[2] else NA
    })

    as.character(na.omit(res))
  }

  if (
    is.null(monitor) ||
      (!is.null(monitor_mode) && monitor_mode %in% c("interpretable", "all"))
  ) {
    # Extract all parameters
    all_params <- unique(c(
      extract_names("^\\s*beta"),
      extract_names("^\\s*alpha"),
      extract_names("^\\s*lambda"),
      extract_names("^\\s*tau"),
      extract_names("^\\s*rho"),
      extract_names("^\\s*sigma"),
      extract_names("^\\s*z_"),
      extract_names("(^|\\W)p_"),
      extract_names("(^|\\W)psi"),
      extract_names("^\\s*r_"),
      extract_names("^\\s*cutpoint")
    ))

    # Remove tau_obs_* (deterministic constants, not stochastic parameters)
    all_params <- all_params[!grepl("^tau_obs", all_params)]

    if (!is.null(monitor_mode) && monitor_mode == "interpretable") {
      # Filter to interpretable parameters only
      monitor <- all_params[
        (grepl("^alpha", all_params) &
          gsub("^alpha_?", "", all_params) %in% response_vars &
          !gsub("^alpha_?", "", all_params) %in% mag_exogenous_vars) |
          grepl("^beta", all_params) |
          grepl("^rho", all_params) |
          grepl("^sigma", all_params) |
          grepl("^psi", all_params) |
          grepl("^p_", all_params) |
          grepl("^z_", all_params) |
          grepl("^r_", all_params) |
          grepl("^cutpoint", all_params) |
          grepl("^K_", all_params) |
          (grepl("^lambda", all_params) &
            gsub("^lambda_?", "", all_params) %in% response_vars &
            !gsub("^lambda_?", "", all_params) %in% mag_exogenous_vars)
      ]
    } else {
      # monitor_mode == "all" or NULL: include everything
      monitor <- all_params

      # Also include response variables (for imputation inspection)
      response_vars_all <- unique(sapply(equations, function(eq) {
        all.vars(formula(eq)[[2]])
      }))

      if (length(response_vars_all) > 0) {
        adj_response_vars <- unlist(lapply(response_vars_all, function(v) {
          get_monitor_vars_hook(family_obj, v)
        }))
        monitor <- unique(c(monitor, adj_response_vars))
      }
    }
  }

  # Add custom monitors provided by the user
  if (length(custom_monitors) > 0) {
    if (!is.null(monitor)) {
      monitor <- unique(c(monitor, custom_monitors))
    } else {
      monitor <- custom_monitors
    }
  }

  # Guarantee monitor is not an empty character vector
  if (is.null(monitor) || length(monitor) == 0) {
    if (exists("all_params") && length(all_params) > 0) {
      monitor <- all_params
    }
    if (is.null(monitor) || length(monitor) == 0) {
      if (length(lines) > 0) {
        assign_lines <- grep("(<-|~)", lines, value = TRUE)
        matches <- regmatches(
          assign_lines,
          regexec("^\\s*([a-zA-Z0-9_.]+)", assign_lines)
        )
        found <- sapply(matches, function(m) if (length(m) >= 2) m[2] else NA)
        monitor <- unique(as.character(na.omit(found)))
      }
    }
  }

  # --- NIMBLE pre-processing ---
  # Ensure all variables used as precision/covariance matrices are numeric matrices
  if (engine == "nimble" && is.list(data)) {
    prec_vars <- names(data)[grepl("^Prec_", names(data))]
    for (pv in prec_vars) {
      if (!is.matrix(data[[pv]])) {
        data[[pv]] <- as.matrix(data[[pv]])
      }
      # [STABILITY] Add small diagonal jitter (nugget) to ensure positive definiteness
      # and prevent numerical singularities during C++ compilation.
      # Ref: Rasmussen & Williams (2006), Gaussian Processes for Machine Learning.
      diag(data[[pv]]) <- diag(data[[pv]]) + 1e-6
    }
  }

  # Add pointwise log-likelihood monitoring if WAIC requested
  if (WAIC) {
    log_lik_params <- unique(extract_names("^\\s*log_lik"))

    if (length(log_lik_params) > 0) {
      monitor <- unique(c(monitor, log_lik_params))
      if (!quiet) {
        message(
          "Monitoring ",
          length(log_lik_params),
          " pointwise log-likelihood parameter(s) for WAIC"
        )
      }
    }
  }

  # Add response variables
  matches <- regmatches(
    model_string,
    gregexpr(
      "\\b([a-zA-Z0-9_.]+)\\s*\\[1:N\\]\\s*~",
      model_string,
      perl = TRUE
    )
  )[[1]]

  m_response_vars <- unique(gsub("\\s*\\[1:N\\]\\s*~", "", matches))
  for (v in m_response_vars) {
    if (!is.null(variability) && v %in% names(variability_list)) {
      next
    }

    if (!v %in% names(data)) {
      base <- sub("[0-9]+$", "", v)
      if (base %in% names(data)) data[[v]] <- data[[base]]
    }
  }

  # Custom inits (e.g. occupancy latent states)
  extension_inits <- get_inits_hook(family_obj, data)

  # Clean up data list: Remove variables not present in the model code to avoid warnings
  vars_to_remove <- character(0)
  keep_vars <- c("zeros")

  for (v in names(data)) {
    if (v %in% keep_vars) {
      next
    }

    if (!grepl(paste0("\\b", v, "\\b"), model_string, perl = TRUE)) {
      vars_to_remove <- c(vars_to_remove, v)
    }
  }

  if (length(vars_to_remove) > 0) {
    for (v in vars_to_remove) {
      data[[v]] <- NULL
    }
  }

  list(
    monitor         = monitor,
    data            = data,
    extension_inits = extension_inits
  )
}
