#' Internal helpers for JAGS/BUGS model generation
#'
#' Utilities for node deduplication, prior retrieval, hierarchical level
#' inference, index calculation, and loop bounding.
#'
#' @keywords internal
NULL

#' Helper: returns b if a is NULL or if a is a list element that doesn't exist
#' @noRd
`%||%` <- function(a, b) {
  tryCatch(if (!is.null(a)) a else b, error = function(e) b)
}

#' Add lines to JAGS model while preventing node redefinition
#'
#' @param current_lines Character vector of current model lines
#' @param new_lines Character vector of new lines to add
#' @param declared_nodes Optional character vector or environment tracking declared node names
#' @return Updated character vector of model lines
#' @keywords internal
safe_add_lines <- function(current_lines, new_lines, declared_nodes = NULL) {
  if (length(new_lines) == 0) return(current_lines)
  
  new_lines <- unlist(strsplit(new_lines, "\n"))
  clean_lines <- c()
  
  is_env <- is.environment(declared_nodes)
  known_nodes <- if (is_env) declared_nodes$nodes else declared_nodes
  if (is.null(known_nodes)) known_nodes <- character(0)
  
  for (line in new_lines) {
    trimmed <- trimws(line)
    if (trimmed == "" || startsWith(trimmed, "#")) {
      clean_lines <- c(clean_lines, line)
      next
    }
    
    if (grepl("[~]|<-", trimmed)) {
      parts <- strsplit(trimmed, "(\\~|<-)")[[1]]
      if (length(parts) > 1) {
        target <- trimws(parts[1])
        target <- sub("\\[.*\\]", "", target)
        
        if (target %in% known_nodes) {
          next
        }
        known_nodes <- c(known_nodes, target)
        if (is_env) declared_nodes$nodes <- known_nodes
      }
    }
    clean_lines <- c(clean_lines, line)
  }
  return(c(current_lines, clean_lines))
}

#' Get prior for a parameter (custom override or default)
#'
#' @param param_name Name of the parameter
#' @param type Type of parameter ("beta", "alpha", etc.)
#' @param default Default prior string
#' @param priors Named list of custom priors
#' @return Formatted JAGS prior statement
#' @keywords internal
get_prior <- function(param_name, type = "beta", default = NULL, priors = NULL) {
  if (!is.null(priors) && param_name %in% names(priors)) {
    return(paste0(param_name, " ~ ", priors[[param_name]]))
  }
  
  if (!is.null(priors) && type %in% names(priors)) {
    return(paste0(param_name, " ~ ", priors[[type]]))
  }
  
  default_prior <- default %||% "dnorm(0, 0.01)"
  
  if (grepl("^[0-9.]+$", trimws(default_prior))) {
    return(paste0(param_name, " <- ", default_prior))
  }
  
  return(paste0(param_name, " ~ ", default_prior))
}

#' Get precision prior by family
#'
#' @param param_name Name of the parameter
#' @param var_name Name of the variable
#' @param priors Named list of custom priors
#' @param family Named vector/list of response families
#' @return Formatted JAGS precision prior statement
#' @keywords internal
get_precision_prior <- function(param_name, var_name, priors = NULL, family = NULL) {
  if (!is.null(priors) && param_name %in% names(priors)) {
    return(paste0(param_name, " ~ ", priors[[param_name]]))
  }
  var_family <- if (!is.null(family) && var_name %in% names(family)) {
    family[[var_name]]
  } else {
    "gaussian"
  }
  fam_obj <- get_family(var_family)
  return(jags_family_precision_prior(fam_obj, param_name))
}

#' Get multiscale level for a variable
#'
#' @param var Variable name
#' @param h_info Hierarchical information list
#' @param equations Optional list of equations
#' @param latent Optional character vector of latent variables
#' @param categorical_vars Optional character vector of categorical variables
#' @return Character level name or NULL
#' @keywords internal
get_var_level <- function(var, h_info, equations = NULL, latent = NULL, categorical_vars = NULL) {
  if (is.null(h_info)) return(NULL)
  
  lvl <- tryCatch(
    infer_variable_level(
      var, 
      h_info$levels, 
      data = h_info$data, 
      equations = equations, 
      latent = latent,
      hierarchy = h_info$hierarchy
    ),
    error = function(e) NULL
  )
  if (!is.null(lvl)) return(lvl)

  lvl_names <- names(h_info$levels)
  for (lvl_name in lvl_names) {
    if (tolower(var) == tolower(lvl_name)) return(lvl_name)
  }

  if (!is.null(categorical_vars)) {
    for (cat_name in names(categorical_vars)) {
      if (
        var %in% categorical_vars[[cat_name]]$dummies ||
        var == paste0(cat_name, "_dummy")
      ) {
        return(get_var_level(cat_name, h_info, equations = equations, latent = latent))
      }
    }
  }

  return(NULL)
}

#' Get hierarchical level name for a random effect grouping variable
#' @keywords internal
get_random_level <- function(response, group_var, h_info, equations = NULL, latent = NULL, categorical_vars = NULL) {
  if (is.null(h_info)) return(NULL)
  return(get_var_level(group_var, h_info, equations = equations, latent = latent, categorical_vars = categorical_vars))
}

#' Get the finest level that is a descendant of all given variables' levels
#' @keywords internal
get_finest_level_of_vars <- function(vars, h_info, equations = NULL, latent = NULL, categorical_vars = NULL) {
  if (is.null(h_info) || is.null(h_info$hierarchy)) return(NULL)

  v_lvls <- unique(na.omit(sapply(vars, function(v) {
    get_var_level(v, h_info, equations = equations, latent = latent, categorical_vars = categorical_vars)
  })))
  if (length(v_lvls) == 0) return(NULL)
  if (length(v_lvls) == 1) return(v_lvls)

  paths <- strsplit(h_info$hierarchy, "\\s*;\\s*")[[1]]
  path_list <- lapply(paths, function(p) {
    trimws(strsplit(p, "\\s*>\\s*")[[1]])
  })
  all_lvls <- unique(unlist(path_list))

  candidates <- c()
  for (lvl in all_lvls) {
    match_all <- TRUE
    for (v_lvl in v_lvls) {
      if (!is_valid_structure_mapping(v_lvl, lvl, h_info, allow_identity = TRUE)) {
        match_all <- FALSE
        break
      }
    }
    if (match_all) {
      candidates <- c(candidates, lvl)
    }
  }

  if (length(candidates) == 0) return(NULL)

  finest <- candidates[1]
  if (length(candidates) > 1) {
    for (i in 2:length(candidates)) {
      if (is_valid_structure_mapping(finest, candidates[i], h_info, allow_identity = TRUE)) {
        finest <- candidates[i]
      }
    }
  }
  return(finest)
}

#' Get level for a structure
#' @keywords internal
get_struct_lvl <- function(s_name, h_info) {
  if (is.null(h_info) || is.null(h_info$structure_levels)) return(NULL)
  return(h_info$structure_levels[[s_name]])
}

#' Get loop bound (N or N_<Level>) for a response variable
#' @keywords internal
get_loop_bound <- function(response, h_info, default_N = "N") {
  if (is.null(h_info)) return("N")
  lvl <- get_var_level(response, h_info)
  if (is.null(lvl)) {
    return(default_N)
  }
  return(paste0("N_", lvl))
}

#' Get zeros vector name (zeros or zeros_<Level>)
#' @keywords internal
get_zeros_name <- function(response, h_info) {
  if (is.null(h_info)) return("zeros")
  lvl <- get_var_level(response, h_info)
  if (is.null(lvl)) {
    if (!is.null(h_info$hierarchy)) {
      h_lvls <- trimws(strsplit(h_info$hierarchy, ">")[[1]])
      return(paste0("zeros_", h_lvls[length(h_lvls)]))
    }
    return("zeros")
  }
  return(paste0("zeros_", lvl))
}

#' Get index expression for accessing a predictor from a coarser level
#' @keywords internal
get_pred_index <- function(pred, response_level, h_info, equations = NULL, latent = NULL, categorical_vars = NULL) {
  if (is.null(h_info)) return(paste0(pred, "[i]"))

  pred_lvl <- get_var_level(pred, h_info, equations = equations, latent = latent, categorical_vars = categorical_vars)
  if (is.null(pred_lvl)) return(paste0(pred, "[i]"))

  if (is.null(response_level) || pred_lvl == response_level) {
    return(paste0(pred, "[i]"))
  }

  idx_name <- paste0(pred_lvl, "_idx_", response_level)
  return(paste0(pred, "[", idx_name, "[i]]"))
}

#' Validate if a structure's level can map to a response's level
#' @keywords internal
is_valid_structure_mapping <- function(s_lvl, r_lvl, h_info, allow_identity = FALSE) {
  if (is.null(r_lvl)) return(TRUE)
  if (is.null(s_lvl)) {
    if (!is.null(h_info) && !is.null(h_info$hierarchy)) return(FALSE)
    return(TRUE)
  }
  if (s_lvl == r_lvl && !allow_identity) return(FALSE)
  if (is.null(h_info) || is.null(h_info$hierarchy)) return(TRUE)

  paths <- strsplit(h_info$hierarchy, "\\s*;\\s*")[[1]]
  for (path in paths) {
    levels <- trimws(strsplit(path, "\\s*>\\s*")[[1]])
    s_idx <- match(s_lvl, levels)
    r_idx <- match(r_lvl, levels)
    if (!is.na(s_idx) && !is.na(r_idx) && s_idx <= r_idx) {
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Get index for mapping structure levels to response levels
#' @keywords internal
get_struct_index <- function(s_name, response, h_info, equations = NULL, latent = NULL, categorical_vars = NULL) {
  s_lvl <- get_struct_lvl(s_name, h_info)
  r_lvl <- get_var_level(response, h_info, equations = equations, latent = latent, categorical_vars = categorical_vars)
  
  if (is.null(s_lvl) || is.null(r_lvl) || s_lvl == r_lvl) {
    return("i")
  }
  if (nchar(s_lvl) == 0 || nchar(r_lvl) == 0) return("i")
  
  return(paste0(s_lvl, "_idx_", r_lvl, "[i]"))
}

#' Check if random effect grouping level is compatible with response level
#' @keywords internal
is_valid_random_level <- function(response, group_var, h_info, family = NULL, equations = NULL, latent = NULL, categorical_vars = NULL, quiet = FALSE) {
  if (is.null(h_info)) return(TRUE)

  resp_lvl <- get_var_level(response, h_info, equations = equations, latent = latent, categorical_vars = categorical_vars)
  grp_lvl <- get_var_level(group_var, h_info, equations = equations, latent = latent, categorical_vars = categorical_vars)

  if (is.null(resp_lvl) || is.null(grp_lvl)) {
    if (is.null(grp_lvl) && !is.null(h_info$data)) {
      for (lvl_name in names(h_info$data)) {
        if (group_var %in% colnames(h_info$data[[lvl_name]])) {
          grp_lvl <- lvl_name
          break
        }
      }
    }
    
    if (is.null(grp_lvl) && !is.null(h_info$hierarchy)) return(FALSE)
    if (is.null(resp_lvl) && !is.null(h_info$hierarchy)) return(FALSE)
    if (is.null(resp_lvl) || is.null(grp_lvl)) return(TRUE)
  }

  dist <- if (!is.null(family) && response %in% names(family)) family[[response]] else "gaussian"
  allow_id <- dist %in% c("poisson", "binomial", "bernoulli")
  
  if (!allow_id && dist == "gaussian" && !is.null(h_info$data) && !is.null(grp_lvl) && !is.null(resp_lvl)) {
    resp_df <- h_info$data[[resp_lvl]]
    if (group_var %in% colnames(resp_df)) {
      n_obs <- nrow(resp_df)
      n_unique <- length(unique(resp_df[[group_var]]))
      if (n_unique < n_obs) {
        allow_id <- TRUE
      } else {
        if (!quiet) {
          warn_id <- paste0("warn_", response, "_", group_var)
          if (!isTRUE(getOption(warn_id))) {
            options(setNames(list(TRUE), warn_id))
            warning(sprintf(
              "Random effect (1|%s) for response '%s' ignored: requires repeated measures for identifiability at the '%s' level.",
              group_var, response, resp_lvl
            ))
          }
        }
        return(FALSE)
      }
    }
  }

  valid <- is_valid_structure_mapping(grp_lvl, resp_lvl, h_info, allow_identity = allow_id)
  
  if (!valid && !quiet && !is.null(h_info$hierarchy)) {
    warn_id <- paste0("warn_invalid_", response, "_", group_var)
    if (!isTRUE(getOption(warn_id))) {
      options(setNames(list(TRUE), warn_id))
      warning(sprintf(
        "Random effect (1|%s) for response '%s' is at an invalid hierarchical level and will be ignored.",
        group_var, response
      ))
    }
  }
  
  return(valid)
}

#' Resolve group indexing for hierarchical random effects
#' @keywords internal
get_group_idx_string <- function(response, r_name, hierarchical_info, default_N = "N") {
  group_idx <- paste0("group_", r_name, "[i]")

  if (!is.null(hierarchical_info)) {
    resp_bound <- get_loop_bound(response, hierarchical_info, default_N = default_N)
    grp_bound <- get_loop_bound(r_name, hierarchical_info, default_N = default_N)

    resp_lvl <- sub("^N_", "", resp_bound)
    grp_lvl <- sub("^N_", "", grp_bound)

    finest_lvl <- trimws(strsplit(hierarchical_info$hierarchy, ">")[[1]])
    finest_lvl <- finest_lvl[length(finest_lvl)]

    if (resp_lvl == "N") resp_lvl <- finest_lvl
    if (grp_lvl == "N") grp_lvl <- finest_lvl

    if (resp_lvl != grp_lvl) {
      group_idx <- paste0(
        "group_",
        r_name,
        "[",
        grp_lvl,
        "_idx_",
        resp_lvl,
        "[i]]"
      )
    }
  }
  return(group_idx)
}

#' Check if a structure is multi-object
#' @keywords internal
is_struct_multi <- function(s_name, hierarchical_info = NULL) {
  if (!is.null(hierarchical_info$structure_multi) && s_name %in% names(hierarchical_info$structure_multi)) {
    status <- hierarchical_info$structure_multi[[s_name]]
    if (!is.null(status)) return(as.logical(status))
  }
  return(FALSE)
}
