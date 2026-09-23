# DAG Layout and Data Assembly
#
# Builds the tidy data frame used by plot_dag() for a single model object.
# Handles: equation extraction, dag_expand_hook, dagitty layout, coordinate
# overrides, node type assignment, edge coefficient matching, and marginal
# effects expansion for multinomial bundles.
#
# Primary function: build_dag_data()
#
# @keywords internal

# Internal helper: null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Build tidy DAG data frame for a single model or formula list
#'
#' @param obj  A `because` model object or a list of formulas.
#' @param label Character label for this model (used in faceting).
#' @param latent Optional character vector of latent variable names.
#' @param family Optional named character vector of family distributions.
#' @param layout Layout algorithm (passed to ggdag::tidy_dagitty).
#' @param coords Optional named list of node coordinates.
#' @param node_size Base node size.
#' @param text_size Node label text size.
#' @param edge_label_size Edge coefficient label size.
#' @param edge_color_scheme One of "directional", "binary", "monochrome".
#' @param type "raw" or "marginal" effects.
#' @param multinomial_probabilities Logical; expand multinomial categories.
#' @param show_coefficients Logical; show edge coefficients.
#' @return A data frame suitable for ggplot2 rendering of a DAG.
#' @keywords internal
build_dag_data <- function(
    obj, label,
    latent = NULL, family = NULL,
    layout = "kk", coords = NULL,
    node_size = 12, text_size = 3.5,
    edge_label_size = 3,
    edge_color_scheme = "directional",
    type = "raw",
    multinomial_probabilities = TRUE,
    show_coefficients = TRUE
) {
    # 1. Extract equations, latent, family, poly_terms from obj
    current_latent    <- latent
    current_family    <- family
    current_poly_terms <- NULL
    equations         <- NULL

    if (inherits(obj, "because")) {
        if (is.null(obj)) return(NULL)
        equations <- if (!is.null(obj$equations)) obj$equations else {
            if (!is.null(obj$input$equations)) obj$input$equations else obj$parameter_map$equations
        }
        if (is.null(equations)) {
            stop("Could not find equations in 'because' object. Please refit the model or manually pass equations.")
        }
        if (is.null(current_latent)) {
            current_latent <- if (!is.null(obj$latent)) obj$latent else obj$input$latent
        }
        if (is.null(current_family)) {
            current_family <- if (!is.null(obj$family)) obj$family else obj$input$family
        }
        current_poly_terms <- if (!is.null(obj$poly_terms)) obj$poly_terms else obj$input$poly_terms
    } else {
        # List of formulas: auto-detect latent variables
        equations <- obj
        if (is.null(current_latent)) {
            lhs_vars <- vapply(equations, function(f) deparse(f[[2]]), character(1))
            rhs_vars <- unique(unlist(lapply(equations, function(f) all.vars(f[[3]]))))
            rhs_only <- setdiff(rhs_vars, lhs_vars)
            latent_pattern <- "^[Ll]([0-9]+|at(ent?)?[0-9]*(_\\w*)?)$"
            auto_latent <- rhs_only[grepl(latent_pattern, rhs_only, perl = TRUE)]
            if (length(auto_latent) > 0) {
                current_latent <- unique(c(current_latent, auto_latent))
            }
        }
    }

    current_equations <- equations

    # 2. S3 Extension Hook: expand compound latent nodes (e.g. occupancy psi_/p_)
    family_obj <- current_family
    if (!is.null(current_family)) {
        for (f in unique(tolower(unlist(current_family)))) {
            class(family_obj) <- unique(c(paste0("because_family_", f), class(family_obj)))
        }
    }
    dag_expansion    <- dag_expand_hook(family_obj, current_equations, current_latent)
    current_equations <- dag_expansion$equations
    current_latent    <- dag_expansion$latent
    compound_groups   <- dag_expansion$compound_groups
    extra_dag_edges   <- dag_expansion$extra_edges %||% character()

    # 3. Convert equations to dagitty string (or use cached dag from fitted model)
    induced_cors <- NULL
    if (inherits(obj, "because")) {
        induced_cors <- if (!is.null(obj$induced_correlations)) obj$induced_correlations else obj$input$induced_correlations
    }

    dag_obj          <- NULL
    interaction_nodes <- NULL

    if (inherits(obj, "because") && !is.null(obj$dag)) {
        dag_obj <- obj$dag
    } else {
        dag_result <- equations_to_dag_string(
            current_equations,
            induced_cors,
            family           = current_family,
            poly_terms       = current_poly_terms,
            collapse_expanded = (type == "marginal"),
            extra_edges      = extra_dag_edges
        )
        dag_str           <- dag_result$dag_string
        interaction_nodes <- dag_result$interaction_nodes
        dag_obj           <- dagitty::dagitty(dag_str)
    }

    # 4. Layout: tidy_dagitty (with optional multiscale tier layout)
    if (layout == "multiscale" && inherits(obj, "because") && !is.null(obj$hierarchical_info)) {
        h_info    <- obj$hierarchical_info
        h_paths   <- strsplit(h_info$hierarchy, "\\s*;\\s*")[[1]]
        all_levels <- unique(trimws(unlist(strsplit(h_paths, "\\s*>\\s*"))))

        var_to_lvl <- list()
        for (lvl in names(h_info$levels)) {
            for (v in h_info$levels[[lvl]]) var_to_lvl[[v]] <- lvl
        }

        all_nodes_in_dag <- names(dag_obj)
        new_coords       <- list()
        lvl_depths       <- vapply(all_levels, function(l) get_level_depth(l, h_info$hierarchy), numeric(1))

        for (lvl in all_levels) {
            lvl_nodes <- all_nodes_in_dag[vapply(all_nodes_in_dag, function(n) {
                clean_n <- if (!is.null(interaction_nodes) && n %in% names(interaction_nodes)) {
                    pvars <- all.vars(parse(text = n))
                    if (length(pvars) > 0) pvars[[1]] else n
                } else n
                (var_to_lvl[[clean_n]] %||% "unknown") == lvl
            }, logical(1))]

            if (length(lvl_nodes) > 0) {
                y_val  <- -lvl_depths[lvl]
                x_vals <- seq(-1, 1, length.out = length(lvl_nodes))
                if (length(lvl_nodes) == 1) x_vals <- 0
                for (j in seq_along(lvl_nodes)) {
                    new_coords[[lvl_nodes[j]]] <- c(x_vals[j], y_val)
                }
            }
        }

        unassigned <- setdiff(all_nodes_in_dag, names(new_coords))
        if (length(unassigned) > 0) {
            y_bottom <- min(sapply(new_coords, function(v) v[2])) - 1
            x_vals   <- seq(-1, 1, length.out = length(unassigned))
            for (j in seq_along(unassigned)) {
                new_coords[[unassigned[j]]] <- c(x_vals[j], y_bottom)
            }
        }

        if (is.null(coords)) coords <- list()
        for (n in names(new_coords)) {
            if (is.null(coords[[n]])) coords[[n]] <- new_coords[[n]]
        }
        tidy_dag <- ggdag::tidy_dagitty(dag_obj, layout = "nicely")
    } else {
        tidy_dag <- ggdag::tidy_dagitty(dag_obj, layout = layout)
    }

    dag_data <- as.data.frame(dplyr::as_tibble(tidy_dag))

    # 5. Apply coordinate overrides
    if (!is.null(coords)) {
        if (!is.list(coords)) {
            stop("coords must be a named list of numeric vectors, e.g. list(node = c(x, y)).")
        }
        for (coord_key in names(coords)) {
            new_pos <- coords[[coord_key]]
            nm <- if (coord_key %in% dag_data$name) {
                coord_key
            } else {
                sanitize_term_name(coord_key)
            }
            is_source <- which(dag_data$name == nm)
            if (length(is_source) > 0) {
                dag_data$x[is_source] <- new_pos[1]
                dag_data$y[is_source] <- new_pos[2]
            }
            if ("to" %in% names(dag_data)) {
                is_target <- which(dag_data$to == nm)
                if (length(is_target) > 0) {
                    dag_data$xend[is_target] <- new_pos[1]
                    dag_data$yend[is_target] <- new_pos[2]
                }
            }
        }
    }

    # 6. Compound group column (for bounding box rendering)
    dag_data$compound_group <- NA_character_
    if (!is.null(compound_groups) && length(compound_groups) > 0) {
        for (grp_name in names(compound_groups)) {
            grp <- compound_groups[[grp_name]]
            dag_data$compound_group[dag_data$name %in% grp$nodes] <- grp_name
        }
    }

    # 7. Node label processing
    dag_data$label_display <- gsub("_", "\n", dag_data$name)
    if (length(interaction_nodes) > 0) {
        for (iname in names(interaction_nodes)) {
            dag_data$label_display[dag_data$name == iname] <- interaction_nodes[[iname]]
        }
    }
    if (!is.null(compound_groups) && length(compound_groups) > 0) {
        for (grp_name in names(compound_groups)) {
            grp <- compound_groups[[grp_name]]
            if (!is.null(grp$labels)) {
                for (node_lbl in names(grp$labels)) {
                    dag_data$label_display[dag_data$name == node_lbl] <- grp$labels[[node_lbl]]
                }
            }
        }
    }

    # 8. Dynamic node size based on longest label
    max_chars        <- max(nchar(unlist(strsplit(dag_data$label_display, "\n"))))
    calc_size        <- min(25, 6 + (max_chars * 1.3 * (text_size / 3.5)))
    current_node_size <- if (node_size == 12) max(10, calc_size) else node_size

    # 9. Initialize edge columns
    dag_data$edge_type   <- NA_character_
    dag_data$weight_abs  <- 1.0
    dag_data$val         <- NA_real_
    dag_data$edge_label  <- NA_character_
    dag_data$significant <- NA
    dag_data$curvature   <- 0

    # 10. Node type assignment (Observed / Latent / Interaction)
    dag_data$type <- "Observed"
    if (!is.null(current_latent)) {
        dag_data$type[dag_data$name %in% current_latent] <- "Latent"
    }
    if (length(interaction_nodes) > 0) {
        dag_data$type[dag_data$name %in% names(interaction_nodes)] <- "Interaction"
    }

    # 11. Edge coefficient matching
    edges_meta <- dagitty::edges(dag_obj)
    stats      <- NULL
    quantiles  <- NULL
    me_table   <- NULL

    if (inherits(obj, "because") && !is.null(obj$summary)) {
        stats     <- obj$summary$statistics
        quantiles <- obj$summary$quantiles
        if (type == "marginal" && !is.null(obj$parameter_map)) {
            me_table <- marginal_effects(obj, samples = 1000, multinomial_probabilities = multinomial_probabilities)
        }
    }

    if (nrow(edges_meta) > 0) {
        edges_to_remove <- c()
        expanded_edges  <- list()

        all_names <- if (inherits(dag_data, "tidy_dagitty")) names(dag_data$data) else names(dag_data)
        edge_rows <- if ("to" %in% all_names) which(!is.na(dag_data$to)) else c()

        for (idx in edge_rows) {
            v <- as.character(dag_data$name[idx])
            w <- as.character(dag_data$to[idx])

            meta_match <- which(
                (as.character(edges_meta$v) == v & as.character(edges_meta$w) == w) |
                (as.character(edges_meta$v) == w & as.character(edges_meta$w) == v)
            )

            if (length(meta_match) > 0) {
                m_idx    <- meta_match[1]
                e_type   <- as.character(edges_meta$e[m_idx])
                dag_data$edge_type[idx] <- e_type

                actual_v <- as.character(edges_meta$v[m_idx])
                actual_w <- as.character(edges_meta$w[m_idx])

                processed_as_marginal <- FALSE
                if (any(type == "marginal") && e_type == "->") {
                    if (!is.null(me_table)) {
                        raw_v <- trimws(gsub("[`]", "", gsub("^(psi_|p_|z_)", "", actual_v)))
                        raw_w <- trimws(gsub("[`]", "", gsub("^(psi_|p_|z_)", "", actual_w)))

                        me_resp_names <- trimws(gsub("[`]", "", as.character(me_table$Response)))
                        me_pred_names <- trimws(gsub("[`]", "", as.character(me_table$Predictor)))

                        clean_v <- raw_v
                        if (!(clean_v %in% me_pred_names)) {
                            clean_v <- gsub("(_L|_Q|_C|_dummy|_\\d+|\\[[0-9]+\\])$", "", raw_v)
                        }
                        clean_w <- raw_w
                        if (!(clean_w %in% me_resp_names)) {
                            clean_w <- gsub("(_L|_Q|_C|_dummy|_\\d+|\\[[0-9]+\\])$", "", raw_w)
                        }

                        me_rows <- me_table[me_resp_names == clean_w & me_pred_names == clean_v, ]

                        if (nrow(me_rows) == 0 && (grepl("Income", clean_v) || grepl("Satsfied", clean_w))) {
                            me_rows <- me_table[
                                grepl(substr(clean_w, 1, 7), me_resp_names, ignore.case = TRUE) &
                                grepl(substr(clean_v, 1, 7), me_pred_names, ignore.case = TRUE), ]
                        }

                        if (nrow(me_rows) > 0) {
                            processed_as_marginal <- TRUE
                            edges_to_remove <- unique(c(edges_to_remove, idx))
                            bundle_size     <- nrow(me_rows)
                            curvatures      <- if (bundle_size == 1) 0 else seq(-0.3, 0.3, length.out = bundle_size)

                            for (k in seq_len(bundle_size)) {
                                new_edge <- dag_data[idx, , drop = FALSE]
                                e_val    <- me_rows$Effect[k]
                                e_low    <- me_rows$Lower[k]
                                e_upp    <- me_rows$Upper[k]
                                e_cat    <- as.character(me_rows$Category[k])

                                new_edge$val        <- e_val
                                new_edge$weight_abs <- abs(e_val)
                                new_edge$edge_label <- if (!is.na(e_cat) && e_cat != "NA") {
                                    paste0(e_cat, ": ", round(e_val, 2))
                                } else round(e_val, 2)
                                new_edge$curvature <- curvatures[k]

                                if (edge_color_scheme != "monochrome") {
                                    if (sign(e_low) == sign(e_upp)) {
                                        new_edge$significant <- if (edge_color_scheme == "directional") (if (e_val > 0) "pos" else "neg") else "sig"
                                    } else new_edge$significant <- "ns"
                                } else new_edge$significant <- "default"

                                expanded_edges[[length(expanded_edges) + 1]] <- new_edge
                            }
                        }
                    }
                }

                if (!processed_as_marginal) {
                    if (!is.null(stats)) {
                        val      <- NA
                        sig_cat  <- "default"
                        pname    <- NULL

                        if (e_type == "->") {
                            try_pname <- paste0("beta_", actual_w, "_", actual_v)
                            if (try_pname %in% rownames(stats)) pname <- try_pname
                        } else if (e_type == "<->") {
                            pname1 <- paste0("rho_", v, "_", w)
                            pname2 <- paste0("rho_", w, "_", v)
                            if (pname1 %in% rownames(stats)) pname <- pname1 else if (pname2 %in% rownames(stats)) pname <- pname2
                        }

                        if (!is.null(pname) && pname %in% rownames(stats)) {
                            val     <- stats[pname, "Mean"]
                            sig_cat <- "default"
                            if (!is.null(quantiles) && pname %in% rownames(quantiles)) {
                                lower <- quantiles[pname, "2.5%"]
                                upper <- quantiles[pname, "97.5%"]
                                if (sign(lower) == sign(upper)) {
                                    sig_cat <- if (edge_color_scheme == "directional") (if (val > 0) "pos" else "neg") else "sig"
                                } else sig_cat <- "ns"
                            }

                            if (!is.na(val)) {
                                dag_data$val[idx]         <- val
                                dag_data$weight_abs[idx]  <- abs(val)
                                dag_data$edge_label[idx]  <- round(val, 2)
                                dag_data$significant[idx] <- sig_cat
                            }
                        }
                    }
                }
            }
        }

        if (length(edges_to_remove) > 0) {
            dag_data <- dag_data[-edges_to_remove, , drop = FALSE]
        }
        if (length(expanded_edges) > 0) {
            dag_data <- rbind(dag_data, do.call(rbind, expanded_edges))
        }
    }

    # 12. Fill defaults
    dag_data$weight_abs[is.na(dag_data$weight_abs) & !is.na(dag_data$to)] <- 1.0
    dag_data$edge_type[is.na(dag_data$edge_type) & !is.na(dag_data$to)]   <- "->"

    if (!"significant" %in% names(dag_data)) dag_data$significant <- "default"
    dag_data$significant[is.na(dag_data$significant)] <- "default"
    dag_data$significant <- factor(
        as.character(dag_data$significant),
        levels = c("pos", "neg", "sig", "ns", "default")
    )

    dag_data$final_node_size <- current_node_size
    dag_data$model_label     <- label

    dag_data
}
