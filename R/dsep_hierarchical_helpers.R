# Hierarchical Helpers for D-Separation Testing
#
# Pure utility functions for multi-level / hierarchical model d-sep:
#   - get_var_level_dsep()             find which hierarchy level a variable belongs to
#   - are_levels_orthogonal()          check if two levels are in separate hierarchy chains
#   - is_cross_hierarchy_test()        detect trivially-satisfied cross-branch tests
#   - get_inherited_random_terms()     resolve random effects from ancestral levels
#   - is_valid_structure_mapping_dsep() check if a structure level maps to a response level
#
# @keywords internal

# Get the level name for a variable
get_var_level_dsep <- function(var, hierarchical_info) {
  if (is.null(hierarchical_info) || is.null(hierarchical_info$levels)) {
    return(NULL)
  }
  for (lvl in names(hierarchical_info$levels)) {
    if (var %in% hierarchical_info$levels[[lvl]]) {
      return(lvl)
    }
  }
  return(NULL)
}

# Determine if two hierarchy levels belong to separate, non-nested branches.
# The hierarchy string has the form "a > b > c; d > e" meaning the semicolon
# separates completely independent chains. Two levels are "orthogonal" if they
# appear in different chains and neither is an ancestor of the other.
are_levels_orthogonal <- function(lvl_a, lvl_b, hierarchy_str) {
  if (is.null(lvl_a) || is.null(lvl_b) || lvl_a == lvl_b) return(FALSE)
  paths <- strsplit(hierarchy_str, "\\s*;\\s*")[[1]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1]]))
  # Find which paths each level belongs to
  path_a <- which(sapply(paths, function(p) lvl_a %in% p))
  path_b <- which(sapply(paths, function(p) lvl_b %in% p))
  # Orthogonal = they are in different chains
  if (length(path_a) == 0 || length(path_b) == 0) return(FALSE)
  !any(path_a %in% path_b)
}

#' Check if a d-sep test is a cross-hierarchy test (trivially satisfied)
#'
#' A test Response _||_ FocalPredictor | \code{\{...\}} is "cross-hierarchy" when the
#' response and the focal predictor live in orthogonal hierarchical branches
#' (e.g., one is species-level, the other is site-level) AND the conditioning
#' set does not contain a variable from a level that connects the two branches
#' (typically the observation level). Such tests are trivially independent by
#' design of the hierarchical data structure and cannot be run in JAGS because
#' no cross-level index exists between orthogonal branches.
#'
#' @param test_eq A formula with attribute "test_var" naming the focal predictor.
#' @param hierarchical_info List with $levels (named list) and $hierarchy (string).
#' @return TRUE if the test should be skipped; FALSE otherwise.
#' @keywords internal
is_cross_hierarchy_test <- function(test_eq, hierarchical_info) {
  if (is.null(hierarchical_info) ||
      is.null(hierarchical_info$levels) ||
      is.null(hierarchical_info$hierarchy)) {
    return(FALSE)
  }
  test_var <- attr(test_eq, "test_var")
  
  # FALLBACK: If attribute is missing, assume first predictor on RHS is the focal one
  if (is.null(test_var)) {
    rhs_vars <- all.vars(test_eq)[-1]
    if (length(rhs_vars) > 0) test_var <- rhs_vars[1]
  }

  if (is.null(test_var)) return(FALSE)

  # Response is the LHS
  resp <- as.character(test_eq)[2]
  resp <- sub("^psi_", "", resp)

  resp_lvl <- get_var_level_dsep(resp, hierarchical_info)
  pred_lvl <- get_var_level_dsep(test_var, hierarchical_info)

  if (is.null(resp_lvl) || is.null(pred_lvl)) return(FALSE)
  if (resp_lvl == pred_lvl) return(FALSE)

  # Are the two levels in separate hierarchy chains?
  if (!are_levels_orthogonal(resp_lvl, pred_lvl, hierarchical_info$hierarchy)) {
    return(FALSE)
  }

  # Even if orthogonal, check whether any conditioning variable bridges them
  # via the obs (or another shared) level. If the conditioning set contains a
  # variable at the obs level (or any level that appears in BOTH chains), the
  # path could be opened. For safety, only skip if conditioning set has NO
  # bridging variable.
  rhs <- as.character(test_eq)[3]
  cond_terms <- trimws(strsplit(rhs, "\\+")[[1]])
  # Strip random effects (1 | X)
  cond_terms <- cond_terms[!grepl("^\\s*1\\s*\\|", cond_terms)]
  # Strip I(...) interaction terms (they are derived, not level-assigning)
  cond_terms <- cond_terms[!grepl("^\\s*I\\(", cond_terms)]
  cond_terms <- trimws(gsub("\\(.*\\)", "", cond_terms))
  cond_terms <- cond_terms[nchar(cond_terms) > 0]

  # Check if any conditioning variable is at a level that appears in both chains
  hierarchy_str <- hierarchical_info$hierarchy
  paths <- strsplit(hierarchy_str, "\\s*;\\s*")[[1]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1]]))
  path_a <- which(sapply(paths, function(p) resp_lvl %in% p))
  path_b <- which(sapply(paths, function(p) pred_lvl %in% p))

  for (cv in cond_terms) {
    cv_lvl <- get_var_level_dsep(cv, hierarchical_info)
    if (is.null(cv_lvl)) next
    cv_in_a <- any(sapply(paths[path_a], function(p) cv_lvl %in% p))
    cv_in_b <- any(sapply(paths[path_b], function(p) cv_lvl %in% p))

    # A variable bridges only if it appears in BOTH lineages. 
    # Example: 'obs' appears in both 'site > survey > obs' and 'species > obs'.
    # A variable at 'site' level does not bridge 'survey' to 'species'.
    if (cv_in_a && cv_in_b) return(FALSE)
  }

  return(TRUE)
}

#' Get inherited random terms for a variable based on hierarchy
#'
#' @param var Variable name
#' @param hierarchical_info List with $levels, $hierarchy, and $link_vars
#' @return List of random terms (group, type)
#' @keywords internal
get_inherited_random_terms <- function(var, hierarchical_info) {
  if (is.null(hierarchical_info) ||
      is.null(hierarchical_info$levels) ||
      is.null(hierarchical_info$hierarchy) ||
      is.null(hierarchical_info$link_vars)) {
    return(list())
  }

  lvl <- get_var_level_dsep(var, hierarchical_info)
  if (is.null(lvl)) return(list())

  # Find all ancestors of this level in the hierarchy string
  # (Splitting by ; handles separate chains)
  paths <- strsplit(hierarchical_info$hierarchy, "\\s*;\\s*")[[1]]
  paths <- lapply(paths, function(p) trimws(strsplit(p, "\\s*>\\s*")[[1]]))

  ancestors <- character(0)
  for (path in paths) {
    idx <- match(lvl, path)
    if (!is.na(idx) && idx > 1) {
      ancestors <- c(ancestors, path[1:(idx - 1)])
    }
  }
  ancestors <- unique(ancestors)

  inherited <- list()
  for (anc in ancestors) {
    if (anc %in% names(hierarchical_info$link_vars)) {
      group_var <- hierarchical_info$link_vars[[anc]]
      if (!is.null(group_var)) {
        inherited[[length(inherited) + 1]] <- list(
          group = group_var,
          type = "intercept"
        )
      }
    }
  }

  return(inherited)
}

# Check if a structure's level can map to a response's level
is_valid_structure_mapping_dsep <- function(s_lvl, r_lvl, hierarchical_info) {
  if (is.null(s_lvl) || is.null(r_lvl) || s_lvl == r_lvl) {
    return(TRUE)
  }
  if (is.null(hierarchical_info) || is.null(hierarchical_info$hierarchy)) {
    return(TRUE)
  }

  paths <- strsplit(hierarchical_info$hierarchy, "\\s*;\\s*")[[1]]
  for (path in paths) {
    levels <- trimws(strsplit(path, "\\s*>\\s*")[[1]])
    s_idx <- match(s_lvl, levels)
    r_idx <- match(r_lvl, levels)
    # Structure must be at or above response in the hierarchy
    if (!is.na(s_idx) && !is.na(r_idx) && s_idx <= r_idx) {
      return(TRUE)
    }
  }
  return(FALSE)
}
