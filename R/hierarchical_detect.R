#' Detect, Validate, and Prepare Hierarchical and Formula Components
#'
#' @noRd
detect_and_validate_hierarchical <- function(
  data,
  equations,
  fixed_eqs_temp,
  parsed_random_temp,
  global_random_vars,
  levels,
  multiscale,
  hierarchy,
  link_vars,
  latent,
  random,
  structure,
  structure_multi = NULL,
  structure_levels = NULL,
  id_col = NULL,
  quiet = FALSE
) {
  # Data is hierarchical if it's a list (not dataframe)
  is_list_data <- is.list(data) && !is.data.frame(data)
  is_hierarchical <- FALSE
  hierarchical_info <- NULL

  if (is_list_data) {
    # Get all variables from fixed equations for auto-detection
    eq_vars <- unique(unlist(lapply(fixed_eqs_temp, all.vars)))

    # Exclude variables used in random terms or global random arg
    all_random_groups <- unique(c(
      vapply(parsed_random_temp$random_terms, function(x) x$group, character(1)),
      global_random_vars
    ))
    eq_vars <- setdiff(eq_vars, all_random_groups)

    # Auto-detect if levels not provided
    if (is.null(levels)) {
      auto_result <- auto_detect_hierarchical(data, eq_vars, quiet = quiet)
      levels <- auto_result$levels

      if (is.null(multiscale)) {
        multiscale <- auto_result$hierarchy
      }
      if (is.null(link_vars)) {
        link_vars <- auto_result$link_vars
      }
    }

    # Validate with provided or auto-detected values
    if (!is.null(levels) && length(levels) > 0) {
      is_hierarchical <- TRUE

      if (is.null(multiscale) && !is.null(hierarchy)) multiscale <- hierarchy

      # Identify deterministic responses (LHS of I() equations)
      det_responses <- character(0)
      for (eq in equations) {
        raw_eq <- formula(eq)
        if (length(raw_eq) >= 3) {
          resp <- as.character(raw_eq[[2]])
          rhs_str <- paste(deparse(raw_eq[[3]]), collapse = " ")
          if (grepl("^I\\(", trimws(rhs_str))) det_responses <- c(det_responses, resp)
        }
      }

      validate_hierarchical_data(
        data,
        levels,
        multiscale,
        link_vars,
        latent_vars = latent,
        equations = equations,
        deterministic_vars = det_responses
      )

      # Infer hierarchy from random effects if still not set
      if (is.null(multiscale)) {
        multiscale <- parse_hierarchy_from_random(random, data)
        if (is.null(multiscale)) {
          stop(
            "Multiscale data detected but 'multiscale' not specified. ",
            "Provide either:\n",
            "  1. 'multiscale' argument (e.g., \"site_year > individual\"), or\n",
            "  2. Nested random effects (e.g., ~(1|site/individual))"
          )
        }
      }

      hierarchical_info <- list(
        data = data,
        levels = levels,
        hierarchy = multiscale,
        link_vars = link_vars,
        latent_vars = latent,
        deterministic_vars = det_responses
      )

      if (!is.null(structure_multi)) {
        hierarchical_info$structure_multi <- structure_multi
      }
      if (!is.null(structure_levels)) {
        hierarchical_info$structure_levels <- structure_levels
      }

      if (!quiet) {
        message("Multiscale causal structure detected: ", multiscale)
      }
    }
  } else {
    if (!is.null(levels) && (!is.null(multiscale) || !is.null(hierarchy))) {
      is_hierarchical <- TRUE
      hierarchical_info <- list(
        data = data,
        levels = levels,
        hierarchy = if (!is.null(multiscale)) multiscale else hierarchy,
        link_vars = link_vars
      )

      if (!is.null(structure_multi)) {
        hierarchical_info$structure_multi <- structure_multi
      }
      if (!is.null(structure_levels)) {
        hierarchical_info$structure_levels <- structure_levels
      }
    }
  }

  row_ids <- NULL
  if (is_hierarchical && !is.null(id_col) && is.list(data) && !is.data.frame(data)) {
    for (.lvl in names(data)) {
      .lvl_df <- data[[.lvl]]
      if (is.data.frame(.lvl_df) && id_col %in% names(.lvl_df)) {
        .col <- .lvl_df[[id_col]]
        if (is.character(.col) || is.factor(.col)) {
          row_ids <- as.character(.col)
          break
        }
      }
    }
  }

  if (is_hierarchical && !is.null(structure) && is.null(hierarchical_info$structure_levels)) {
    hierarchical_info$structure_levels <- auto_detect_structure_levels(structure, hierarchical_info, quiet = quiet)
  }

  list(
    is_hierarchical   = is_hierarchical,
    hierarchical_info = hierarchical_info,
    levels            = levels,
    multiscale        = multiscale,
    link_vars         = link_vars,
    row_ids           = row_ids
  )
}


#' Parse Model Random Effects from Equations and Deprecated Argument
#'
#' @noRd
parse_model_random_effects <- function(equations, random = NULL) {
  parsed_random <- extract_random_effects(equations)
  equations <- parsed_random$fixed_equations
  random_terms <- parsed_random$random_terms

  if (!is.null(random)) {
    resp_vars <- vapply(
      parsed_random$fixed_equations,
      function(eq) as.character(eq[[2]]),
      character(1)
    )
    rand_char  <- as.character(random)
    rand_rhs   <- rand_char[length(rand_char)]
    rand_parts <- trimws(strsplit(rand_rhs, "\\+")[[1]])
    rand_parts <- rand_parts[grepl("\\|", rand_parts)]
    rand_inline <- paste(rand_parts, collapse = " + ")
    migration_lines <- vapply(
      resp_vars,
      function(r) paste0("  ", r, " ~ ... + ", rand_inline),
      character(1)
    )
    warning(
      "The 'random' argument is deprecated and will be removed in a future ",
      "version. Embed random effects directly in the equation formulas instead:\n",
      paste(migration_lines, collapse = "\n"), "\n",
      "See vignette('08_multiscale_models') for details.",
      call. = FALSE
    )

    global_random_terms <- parse_global_random(random, equations)
    random_terms <- c(random_terms, global_random_terms)

    if (length(random_terms) > 0) {
      keys <- vapply(
        random_terms,
        function(x) paste(x$response, x$group, sep = "|"),
        character(1)
      )
      random_terms <- random_terms[!duplicated(keys)]
    }
  }

  list(
    equations    = equations,
    random_terms = random_terms
  )
}


#' Process and Expand Polynomial Terms in Formulas
#'
#' @noRd
process_polynomial_terms <- function(
  equations,
  is_hierarchical = FALSE,
  levels = NULL,
  hierarchical_info = NULL,
  quiet = FALSE
) {
  all_poly_terms <- get_all_polynomial_terms(equations)

  if (!is.null(all_poly_terms)) {
    equations <- lapply(equations, function(eq) {
      poly_terms <- extract_polynomial_terms(eq)
      expand_polynomial_formula(eq, poly_terms)
    })

    if (!quiet) {
      message(
        "Detected ",
        length(all_poly_terms),
        " polynomial term(s): ",
        paste(
          vapply(all_poly_terms, function(x) x$original, character(1)),
          collapse = ", "
        )
      )
    }

    if (is_hierarchical && !is.null(levels)) {
      for (poly in all_poly_terms) {
        base_var <- poly$base_var
        new_var <- poly$internal_name

        for (lvl_name in names(levels)) {
          if (base_var %in% levels[[lvl_name]]) {
            levels[[lvl_name]] <- c(levels[[lvl_name]], new_var)
            break
          }
        }
      }

      if (!is.null(hierarchical_info)) {
        hierarchical_info$levels <- levels
      }
    }
  }

  list(
    equations         = equations,
    all_poly_terms    = all_poly_terms,
    levels            = levels,
    hierarchical_info = hierarchical_info
  )
}


#' Assemble Hierarchical Dataset for Model Fitting
#'
#' @noRd
assemble_hierarchical_dataset <- function(
  data,
  hierarchical_info,
  equations,
  random_terms,
  latent,
  all_poly_terms,
  original_raw_data_before_preprocess,
  engine = "jags",
  quiet = FALSE
) {
  eq_vars <- unique(unlist(lapply(equations, all.vars)))

  if (length(random_terms) > 0) {
    random_vars <- unique(vapply(
      random_terms,
      function(x) x$group,
      character(1)
    ))
    eq_vars <- unique(c(eq_vars, random_vars))
  }

  if (!is.null(hierarchical_info) && !is.null(hierarchical_info$link_vars)) {
    eq_vars <- unique(c(eq_vars, unlist(hierarchical_info$link_vars)))
  }

  if (!is.null(attr(data, "categorical_vars"))) {
    cat_vars <- attr(data, "categorical_vars")
    for (cv in names(cat_vars)) {
      if (cv %in% eq_vars) {
        eq_vars <- unique(c(eq_vars, cat_vars[[cv]]$dummies))
      }
    }
  }

  if (!is.null(latent)) {
    eq_vars <- setdiff(eq_vars, latent)
  }

  if (!is.null(all_poly_terms)) {
    base_poly_vars <- vapply(
      all_poly_terms,
      function(x) x$base_var,
      character(1)
    )
    derived_poly_vars <- vapply(
      all_poly_terms,
      function(x) x$internal_name,
      character(1)
    )

    eq_vars <- unique(c(eq_vars, base_poly_vars))
    eq_vars <- setdiff(eq_vars, derived_poly_vars)
  }

  rhs_vars <- unique(unlist(lapply(equations, function(eq) {
    if (length(eq) == 3) {
      all.vars(eq[[3]])
    } else {
      character(0)
    }
  })))

  original_raw_data <- original_raw_data_before_preprocess

  if (!is.null(hierarchical_info)) {
    hierarchical_info$data <- preprocess_categorical_vars(
      hierarchical_info$data,
      dummy_vars = rhs_vars,
      quiet = quiet
    )

    if (is.list(data) && !is.data.frame(data)) {
      data <- hierarchical_info$data
    }
  }

  cat_vars <- attr(data, "categorical_vars")
  if (!is.null(cat_vars)) {
    for (cv_name in names(cat_vars)) {
      if (cv_name %in% rhs_vars) {
        eq_vars <- c(eq_vars, cat_vars[[cv_name]]$dummies)
      }
    }
    eq_vars <- unique(eq_vars)
  }

  prep_res <- prepare_hierarchical_jags_data(hierarchical_info, eq_vars)
  data <- prep_res$data_list

  if (!is.null(cat_vars)) {
    attr(data, "categorical_vars") <- cat_vars
  }

  data <- c(data, prep_res$n_vec)

  if (!quiet) {
    message(
      paste0("Prepared hierarchical data for ", toupper(engine), ": "),
      paste(
        names(prep_res$n_vec),
        unlist(prep_res$n_vec),
        sep = "=",
        collapse = ", "
      )
    )
  }

  original_data <- original_raw_data

  list(
    data              = data,
    hierarchical_info = hierarchical_info,
    original_data     = original_data
  )
}


#' Validate Binomial Response Variables
#'
#' @noRd
validate_binomial_responses <- function(data, equations, family) {
  if (is.null(family)) return(invisible(TRUE))

  responses <- vapply(
    equations,
    function(eq) as.character(formula(eq)[2]),
    character(1)
  )

  for (i in seq_along(responses)) {
    resp <- responses[i]
    fam <- NULL

    if (!is.null(names(family))) {
      if (resp %in% names(family)) {
        fam <- family[[resp]]
      }
    } else {
      if (length(family) == 1) {
        fam <- family
      } else if (length(family) == length(responses)) {
        if (is.list(family)) fam <- family[[i]] else fam <- family[i]
      }
    }

    if (!is.null(fam) && fam == "binomial" && resp %in% names(data)) {
      vals <- data[[resp]]
      if (any(!vals %in% c(0, 1, NA))) {
        unique_vals <- unique(vals[!is.na(vals)])
        stop(sprintf(
          "Binomial response variable '%s' contains invalid values: {%s}. Binomial variables must be strictly 0 or 1.",
          resp,
          paste(head(unique_vals, 5), collapse = ", ")
        ))
      }
    }
  }
  invisible(TRUE)
}
