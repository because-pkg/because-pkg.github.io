#' Plot DAG from Equations or Fitted Model
#'
#' Visualizes the Directed Acyclic Graph (DAG) implied by a set of equations
#' or a fitted `because` model. If a fitted model is provided:
#' * Path coefficients are displayed on the edges.
#' * Edge colors are determined by the `edge_color_scheme` argument (e.g., Red/Blue for directional, Black/Grey for binary).
#' * Structural edges (without fitted data) are colored black.
#' * Edge thickness scales with the absolute effect size.
#'
#' @param x A list of formulas (equations), a `because` model object, or a list of these.
#' @param layout The layout algorithm to use for positioning nodes (default \code{"kk"}, Kamada-Kawai).
#'   Any layout name accepted by \code{\link[ggraph]{create_layout}} can be used, including:
#'   \itemize{
#'     \item \code{"kk"} — Kamada-Kawai spring layout (default; good general-purpose layout).
#'     \item \code{"sugiyama"} — hierarchical/layered layout; useful for simple chains.
#'     \item \code{"fr"} — Fruchterman-Reingold spring layout.
#'     \item \code{"nicely"} — automatic choice based on graph properties.
#'     \item \code{"circle"} — nodes arranged in a circle.
#'   }
#'   Override node positions entirely with the \code{coords} argument.
#' @param latent Character vector of latent variable names. Overrides the model's
#'   latent variables if provided. If \code{NULL} (default) and plotting from formulas, 
#'   latent variables are **automatically detected** if they follow the SEM 
#'   naming convention (e.g., \code{L1}, \code{L2}, \code{Latent1}, \code{lat_foo}) 
#'   and only appear on the RHS of equations.
#' @param node_size Size of the nodes (default 14).
#' @param node_color Color of the node border (default "black").
#' @param node_fill Color of the node interior (default "white").
#' @param node_stroke Thickness of the node border (default 1.5).
#' @param text_size Size of the labels (default 4).
#' @param edge_width_range Vector of length 2 defining the range of arrow widths (min, max) based on effect size.
#' @param edge_color_scheme Character; one of "directional" (default), "binary", or "monochrome".
#' "directional" colors edges red/blue/grey based on effect direction and whether the 95\% CI excludes zero.
#' "binary" colors edges black/grey based on whether the 95\% CI excludes zero (black) or includes it (grey).
#' "monochrome" colors all edges black.
#' @param show_coefficients Logical; whether to print coefficient values on edges (only for fitted models).
#' @param coords Optional named list of coordinates for the nodes, e.g. \code{list(A = c(1, 1), B = c(2, 2))}.
#' If provided, these will override the \code{layout} algorithm. **Partial coordinates** are supported: 
#' nodes not included in \code{coords} will be positioned according to the automatic \code{layout}. 
#' For deterministic nodes (interactions and powers), you can use the original formula string 
#' as the key (e.g., \code{"I(age^2)" = c(x, y)} or \code{"X:Y" = c(x, y)}).
#' @param family Optional named character vector of families for response variables.
#'
#' @return A `ggplot` object that can be further customized with standard ggplot2 functions (e.g., `+ theme_...()`, `+ ggtitle(...)`).
#'
#' @details
#' \strong{Interaction and Deterministic Nodes:}
#' Interaction terms (e.g. \code{BM:M}) and \code{I()} transformations are
#' rendered as **explicit intermediate nodes** (grey diamonds), following the
#' Interaction DAG (IDAG) convention of Attia, Holliday & Oldmeadow (2022).
#' This makes the deterministic nature of these terms visually clear and is
#' consistent with how \code{because_dsep} treats them for d-separation.
#'
#' \strong{Latent Variable Auto-detection:}
#' When plotting from a list of formulas, \code{plot_dag} automatically identifies 
#' latent variables (rendering them as circles) if they match common SEM naming 
#' conventions (like \code{L1}, \code{Latent}, \code{lat_climate}) and never 
#' appear as the response (LHS) of an equation.
#'
#' \strong{Manual Positioning with Formula Strings:}
#' When using \code{coords}, you can specify positions for interaction or 
#' power nodes by using their formula representation as the list key. For example:
#' \code{coords = list(weight = c(0,0), "I(age^2)" = c(1,1), "sex:age" = c(2,2))}. 
#' Any node not specified in the list will maintain its position from the 
#' automatic layout.
#'
#' \strong{Random Effects:}
#' Formula terms containing random effects (e.g., \code{(1|year)}) are 
#' automatically filtered out for the structural DAG visualization.
#'
#' @references
#' Attia, J., Holliday, E., & Oldmeadow, C. (2022). A proposal for capturing
#' interaction and effect modification using DAGs.
#' \emph{International Journal of Epidemiology}, 51(4), 1047--1053.
#' \doi{10.1093/ije/dyac105}
#'
#' @export
#' @importFrom stats terms formula
#'
#' @examples
#' \dontrun{
#' # Basic plotting
#' eq <- list(y ~ x + z, x ~ z)
#' plot_dag(eq)
#'
#' # Custom Layout
#' my_coords <- list(
#'   y = c(1, 1),
#'   x = c(0, 0),
#'   z = c(2, 0)
#' )
#' plot_dag(eq, coords = my_coords)
#' }
#'
plot_dag <- function(
    x,
    layout = "kk",
    latent = NULL,
    node_size = 14,
    node_color = "black",
    node_fill = "white",
    node_stroke = 1.5,
    text_size = 4,
    edge_width_range = c(0.5, 2),
    edge_color_scheme = c("directional", "binary", "monochrome"),
    show_coefficients = TRUE,
    coords = NULL,
    family = NULL,
    type = c("raw", "marginal"),
    multinomial_probabilities = TRUE
) {
    edge_color_scheme <- match.arg(edge_color_scheme)
    type <- match.arg(type)

    # Check dependencies
    if (
        !requireNamespace("dagitty", quietly = TRUE) ||
            !requireNamespace("ggdag", quietly = TRUE) ||
            !requireNamespace("ggraph", quietly = TRUE) ||
            !requireNamespace("ggplot2", quietly = TRUE) ||
            !requireNamespace("dplyr", quietly = TRUE)
    ) {
        stop(
            "Packages 'dagitty', 'ggdag', 'ggraph', 'ggplot2', and 'dplyr' are required for plot_dag."
        )
    }

    # Normalize input to a list of objects
    if (
        inherits(x, "because") ||
            inherits(x, "list") && all(sapply(x, inherits, "formula"))
    ) {
        x <- list(Model = x)
    }

    # Build Tidy DAG Data Frame
    combined_dag_data <- NULL

    for (i in seq_along(x)) {
        obj <- x[[i]]
        label <- names(x)[i]
        if (is.null(label)) {
            label <- paste("Model", i)
        }

        # 1. Extract Equations and Latent Info
        current_latent <- latent
        current_family <- family
        current_poly_terms <- NULL

        if (inherits(obj, "because")) {
            eqs <- obj$equations %||% 
                obj$parameter_map$equations %||%
                obj$input$equations %||%
                obj$parameter_map$equations # Fallback
            if (is.null(eqs)) {
                stop(
                    "Could not find equations in 'because' object. Please refit the model or manually pass equations."
                )
            }
            # Use object's latent vars if not overridden
            if (is.null(current_latent)) {
                current_latent <- obj$input$latent
            }
            # Also get family if not provided
            if (is.null(current_family)) {
                current_family <- obj$input$family
            }
            # Read stored poly_terms so diamond nodes are reconstructed correctly
            # (because() expands I(age^2) -> age_pow2 before storing equations)
            current_poly_terms <- obj$input$poly_terms
        } else {
            # List of formulas
            eqs <- obj

            # Auto-detect latent variables from structural equations:
            # A variable that (a) never appears as a LHS response, and (b) matches
            # common SEM latent naming conventions (L1, L2, Lat, Latent, lat_*, etc.)
            # is automatically treated as latent (rendered as a circle).
            if (is.null(current_latent)) {
                lhs_vars <- vapply(eqs, function(f) deparse(f[[2]]), character(1))
                rhs_vars <- unique(unlist(lapply(eqs, function(f) {
                    all.vars(f[[3]])
                })))
                # Variables only on RHS, never a response
                rhs_only <- setdiff(rhs_vars, lhs_vars)
                # Match common SEM latent naming: L1, L2, Latent, latent1, lat_climate, Lat_foo
                # Matches: L + digits, Lat/lat + optional suffix, Latent/latent + optional digits
                latent_pattern <- "^[Ll]([0-9]+|at(ent?)?[0-9]*(_\\w*)?)$"
                auto_latent <- rhs_only[grepl(latent_pattern, rhs_only, perl = TRUE)]
                if (length(auto_latent) > 0) {
                    current_latent <- unique(c(current_latent, auto_latent))
                }
            }
        }

        # Handle occupancy expansion for visualization
        occ_vars <- c()
        if (!is.null(current_family)) {
            occ_vars <- names(current_family)[current_family == "occupancy"]
            if (length(occ_vars) > 0) {
                # Add p_ and psi_ to latent set so they become circles
                latent_to_add <- unlist(lapply(occ_vars, function(v) {
                    c(paste0("p_", v), paste0("psi_", v))
                }))
                current_latent <- unique(c(current_latent, latent_to_add))
            }
        }

        # 2. Convert to dagitty syntax
        induced_cors <- NULL
        if (inherits(obj, "because")) {
            induced_cors <- obj$induced_correlations %||%
                obj$input$induced_correlations
        }

        dag_result <- equations_to_dag_string(
            eqs,
            induced_cors,
            family = current_family,
            poly_terms = current_poly_terms,
            collapse_expanded = (type == "marginal")
        )
        dag_str <- dag_result$dag_string
        interaction_nodes <- dag_result$interaction_nodes
        dag_obj <- dagitty::dagitty(dag_str)

        # 3. Tidy it up using ggdag
        tidy_dag <- ggdag::tidy_dagitty(dag_obj, layout = layout)

        # Extract data frame for plotting (nodes + edges flattened)
        dag_data <- as.data.frame(dplyr::as_tibble(tidy_dag))

        # Apply coordinates if provided (overwrites layout only for specified nodes)
        if (!is.null(coords)) {
            if (!is.list(coords)) {
                stop(
                    "coords must be a named list of numeric vectors, e.g. list(node = c(x, y))."
                )
            }
        # Apply coordinates if provided (overwrites layout only for specified nodes)
        if (!is.null(coords)) {
            if (!is.list(coords)) {
                stop(
                    "coords must be a named list of numeric vectors, e.g. list(node = c(x, y))."
                )
            }
            # Update positions in the tidy data frame
            for (coord_key in names(coords)) {
                new_pos <- coords[[coord_key]]
                
                # Match the exact key or its sanitized internal version (e.g. "I(age^2)" -> "age_pow2")
                # This allow users to use original formula strings in the coords list.
                nm <- if (coord_key %in% dag_data$name) {
                    coord_key
                } else {
                    # sanitize_term_name matches the logic used to create deterministic nodes
                    sanitize_term_name(coord_key)
                }

                # Update source positions (affecting nodes and start of arrows)
                is_source <- which(dag_data$name == nm)
                if (length(is_source) > 0) {
                    dag_data$x[is_source] <- new_pos[1]
                    dag_data$y[is_source] <- new_pos[2]
                }
                
                # Update target positions (affecting end of arrows)
                if ("to" %in% names(dag_data)) {
                    is_target <- which(dag_data$to == nm)
                    if (length(is_target) > 0) {
                        dag_data$xend[is_target] <- new_pos[1]
                        dag_data$yend[is_target] <- new_pos[2]
                    }
                }
            }
        }
        }

        # Identify occupancy nodes and assign to groups for box drawing
        dag_data$occ_species <- NA_character_
        if (length(occ_vars) > 0) {
            for (v in occ_vars) {
                # Pattern match p_Species and psi_Species
                dag_data$occ_species[grepl(
                    paste0("^(p_|psi_)", v, "$"),
                    dag_data$name
                )] <- v
                # Also include the observation node if it exists
                dag_data$occ_species[dag_data$name == v] <- v
            }
        }

        # Label Processing: Wrap text (replace _ with \n) for regular nodes
        dag_data$label_display <- gsub("_", "\n", dag_data$name)
        # Override for interaction/deterministic nodes: show "BM\u00d7M" etc.
        if (length(interaction_nodes) > 0) {
            for (iname in names(interaction_nodes)) {
                dag_data$label_display[
                    dag_data$name == iname
                ] <- interaction_nodes[[iname]]
            }
        }
        # Override for occupancy latent nodes to just 'p' and 'psi'
        if (length(occ_vars) > 0) {
            dag_data$label_display <- ifelse(
                grepl("^p_", dag_data$name) & !is.na(dag_data$occ_species),
                "p",
                dag_data$label_display
            )
            dag_data$label_display <- ifelse(
                grepl("^psi_", dag_data$name) & !is.na(dag_data$occ_species),
                "psi",
                dag_data$label_display
            )
        }

        # Determine Node Size based on longest label (Heuristic)
        # Assuming circular node, diameter needs to cover the rectangular text block.
        # We find the max characters in any single line of the wrapped labels.
        max_chars <- max(nchar(unlist(strsplit(dag_data$label_display, "\n"))))
        # Heuristic: base size + scaling factor * chars * text_size
        # Typically node_size=14 covers ~3 chars at text_size=4.
        # Reduced multiplier to prevent excessive size in compact/Rmd plots.
        calc_size <- 8 + (max_chars * 1.8 * (text_size / 4))

        # Use simple logic: manual node_size is a baseline, but we ensure it fits.
        # However, user asked to "make the boxes size the same... according to the longest".
        # So we use the calculated max size for ALL nodes stringently.
        # We will take the maximum of the default/user input and the calculated requirement.
        current_node_size <- max(node_size, calc_size)

        # Initialize columns for ALL edges
        dag_data$edge_type <- NA_character_
        dag_data$weight_abs <- 1.0 # Default to 1 (solid/opaque) for structural plots
        dag_data$val <- NA_real_
        dag_data$edge_label <- NA_character_ # Initialize edge_label column
        dag_data$significant <- NA # Logical placeholder
        dag_data$curvature <- 0 # Default curvature (straight)

        # Mark node types (Observed vs Latent vs Interaction)
        dag_data$type <- "Observed"
        if (!is.null(current_latent)) {
            dag_data$type[dag_data$name %in% current_latent] <- "Latent"
        }
        # Interaction/deterministic nodes are distinct from both
        if (length(interaction_nodes) > 0) {
            dag_data$type[
                dag_data$name %in% names(interaction_nodes)
            ] <- "Interaction"
        }

        # 3. Add coefficients if available
        # Get edges metadata from dagitty for finding parameters
        edges_meta <- dagitty::edges(dag_obj) # v, w, e

        stats <- NULL
                quantiles <- NULL
        me_table <- NULL
        if (inherits(obj, "because") && !is.null(obj$summary)) {
            stats <- obj$summary$statistics
            quantiles <- obj$summary$quantiles
            
            if (type == "marginal") {
               # Compute marginal effects (subsampled for speed)
               if (!is.null(obj$parameter_map)) {
                  me_table <- marginal_effects(obj, samples = 500, multinomial_probabilities = multinomial_probabilities)
               }
            }
        }

        # Process edges to find parameters
        if ("to" %in% names(dag_data) && nrow(edges_meta) > 0) {
            edge_rows <- which(!is.na(dag_data$to))
            expanded_edges <- list() # Store multinomial bundles here

            for (idx in edge_rows) {
                v <- dag_data$name[idx]
                w <- dag_data$to[idx]

                # Find corresponding edge via match
                meta_match <- which(
                    (edges_meta$v == v & edges_meta$w == w) |
                        (edges_meta$e == "<->" &
                            edges_meta$v == w &
                            edges_meta$w == v)
                )

                if (length(meta_match) > 0) {
                    m_idx <- meta_match[1]
                    e_type <- edges_meta$e[m_idx]
                    dag_data$edge_type[idx] <- e_type

                    if (!is.null(stats)) {
                        # Default path logic (Standard or First Category)
                        val <- NA
                        sig_cat <- "default"
                        pname <- NULL

                        if (e_type == "->") {
                            # Beta: beta_w_v
                            try_pname <- paste0("beta_", w, "_", v)
                            if (try_pname %in% rownames(stats)) pname <- try_pname
                        } else if (e_type == "<->" ) {
                            pname1 <- paste0("rho_", v, "_", w)
                            pname2 <- paste0("rho_", w, "_", v)
                            if (pname1 %in% rownames(stats)) pname <- pname1 else if (pname2 %in% rownames(stats)) pname <- pname2
                        }

                        # Multinomial / Marginal Effects Check
                        if (type == "marginal" && !is.null(me_table) && e_type == "->") {
                            clean_w <- trimws(gsub("[`]", "", gsub("^(psi_|p_|z_)", "", w)))
                            clean_v <- trimws(gsub("[`]", "", gsub("^(psi_|p_|z_)", "", v)))
                            clean_w <- gsub("(_[A-Za-z0-9]+|\\[[0-9]+\\])$", "", clean_w)
                            clean_v <- gsub("(_[A-Za-z0-9]+|\\[[0-9]+\\])$", "", clean_v)
                            
                            me_rows <- me_table[
                                trimws(gsub("[`]", "", as.character(me_table$Response))) == clean_w & 
                                trimws(gsub("[`]", "", as.character(me_table$Predictor))) == clean_v, 
                            ]

                            if (nrow(me_rows) > 0) {
                                # MULTINOMIAL BUNDLE LOGIC
                                if (nrow(me_rows) > 1) {
                                   # We will create one edge row per category.
                                   # The original row (idx) will be replaced or marked for removal.
                                   # Actually, we keep the original row for the first category and add the rest.
                                   base_row <- dag_data[idx, , drop=FALSE]
                                   
                                   # Calculate curvatures: symmetric bundle around 0
                                   # Formula: seq(-0.3, 0.3, length.out = K)
                                   curvatures <- seq(-0.3, 0.3, length.out = nrow(me_rows))

                                   for (k in seq_len(nrow(me_rows))) {
                                      new_edge <- base_row
                                      e_val <- me_rows$Effect[k]
                                      e_low <- me_rows$Lower[k]
                                      e_upp <- me_rows$Upper[k]
                                      e_cat <- me_rows$Category[k]

                                      new_edge$val <- e_val
                                      new_edge$weight_abs <- abs(e_val)
                                      # Robust Label: Category + Magnitude
                                      new_edge$edge_label <- if (!is.na(e_cat)) paste0(e_cat, ": ", round(e_val, 2)) else round(e_val, 2)
                                      new_edge$curvature <- curvatures[k]
                                      
                                      # Sig color
                                      if (edge_color_scheme != "monochrome") {
                                         if (sign(e_low) == sign(e_upp)) {
                                            new_edge$significant <- if (edge_color_scheme == "directional") (if (e_val > 0) "pos" else "neg") else "sig"
                                         } else {
                                            new_edge$significant <- "ns"
                                         }
                                      } else {
                                         new_edge$significant <- "default"
                                      }
                                      
                                      if (k == 1) {
                                         # Overwrite original row
                                         dag_data[idx, ] <- new_edge
                                      } else {
                                         # Add to expansion list
                                         expanded_edges[[length(expanded_edges) + 1]] <- new_edge
                                      }
                                   }
                                   next # Skip standard processing
                                } else {
                                   # Single category (or default expected score shift)
                                   val <- me_rows$Effect[1]
                                   lower <- me_rows$Lower[1]
                                   upper <- me_rows$Upper[1]
                                   if (edge_color_scheme != "monochrome") {
                                      if (sign(lower) == sign(upper)) {
                                         sig_cat <- if (edge_color_scheme == "directional") (if (val > 0) "pos" else "neg") else "sig"
                                      } else sig_cat <- "ns"
                                   }
                                }
                            }
                        } else if (!is.null(pname)) {
                            val <- stats[pname, "Mean"]
                            if (!is.null(quantiles) && pname %in% rownames(quantiles)) {
                                lower <- quantiles[pname, "2.5%"]
                                upper <- quantiles[pname, "97.5%"]
                                if (sign(lower) == sign(upper)) {
                                    sig_cat <- if (edge_color_scheme == "directional") (if (val > 0) "pos" else "neg") else "sig"
                                } else sig_cat <- "ns"
                            }
                        }

                        if (!is.na(val)) {
                            dag_data$val[idx] <- val
                            dag_data$weight_abs[idx] <- abs(val)
                            dag_data$edge_label[idx] <- round(val, 2)
                            dag_data$significant[idx] <- sig_cat
                        }
                    }
                }
            }
            # Combine Expanded Edges
            if (length(expanded_edges) > 0) {
               dag_data <- rbind(dag_data, do.call(rbind, expanded_edges))
            }
        }

        # Fill defaults for plotting if missing
        dag_data$weight_abs[
            is.na(dag_data$weight_abs) & !is.na(dag_data$to)
        ] <- 1.0
        # Default edge type to -> if not found (robustness)
        dag_data$edge_type[
            is.na(dag_data$edge_type) & !is.na(dag_data$to)
        ] <- "->"
        # Ensure factor levels for significance (Re-apply after adding bundles)
        if (!"significant" %in% names(dag_data)) {
            dag_data$significant <- "default"
        }
        dag_data$significant[is.na(dag_data$significant)] <- "default"
        dag_data$significant <- factor(
            as.character(dag_data$significant),
            levels = c("pos", "neg", "sig", "ns", "default")
        )

        dag_data$final_node_size <- current_node_size
        dag_data$model_label <- label

        if (is.null(combined_dag_data)) {
            combined_dag_data <- dag_data
        } else {
            combined_dag_data <- rbind(combined_dag_data, dag_data)
        }
    }

    # Calculate Bounding Boxes for Occupancy Groups
    occupancy_boxes <- NULL
    if (
        !is.null(combined_dag_data) &&
            any(!is.na(combined_dag_data$occ_species))
    ) {
        # Check if we have coordinates
        if (any(!is.na(combined_dag_data$x))) {
            occupancy_boxes <- combined_dag_data |>
                dplyr::filter(
                    !is.na(occ_species) & grepl("^(p_|psi_)", name)
                ) |>
                dplyr::group_by(model_label, occ_species) |>
                dplyr::summarize(
                    xmin = min(x) - 0.25,
                    xmax = max(x) + 0.25,
                    ymin = min(y) - 0.25,
                    ymax = max(y) + 0.25,
                    .groups = "drop"
                )
        }
    }

    # Ensure type is a factor (Observed, Latent, Interaction)
    combined_dag_data$type <- factor(
        combined_dag_data$type,
        levels = c("Observed", "Latent", "Interaction")
    )

    # Calculate uniform node size for plotting and caps
    # Handle case where combined_dag_data might be missing final_node_size or empty
    if (is.null(combined_dag_data$final_node_size)) {
        uniform_node_size <- node_size
    } else {
        uniform_node_size <- max(
            combined_dag_data$final_node_size,
            na.rm = TRUE
        )
        # Fallback if max is -Inf
        if (!is.finite(uniform_node_size)) uniform_node_size <- node_size
    }

    # Define Caps based on node size
    # Reduced multiplier from 1.5 to 1.15 to decrease white space/margins.
    # This prevents edges from disappearing in compact plots.
    cap_size <- ggraph::circle(uniform_node_size * 1.15, "pt")

    p <- ggplot2::ggplot(
        combined_dag_data,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
    )

    # 0. Shadowed Boxes for Occupancy (Background Layer)
    if (!is.null(occupancy_boxes)) {
        p <- p +
            ggplot2::geom_rect(
                data = occupancy_boxes,
                ggplot2::aes(
                    xmin = xmin,
                    xmax = xmax,
                    ymin = ymin,
                    ymax = ymax
                ),
                inherit.aes = FALSE,
                fill = "grey93",
                color = "grey80",
                alpha = 0.5,
                linewidth = 0.4,
                lty = "dashed"
            ) +
            # Species name label for the box
            ggplot2::geom_text(
                data = occupancy_boxes,
                ggplot2::aes(
                    x = (xmin + xmax) / 2,
                    y = ymax + 0.05,
                    label = occ_species
                ),
                inherit.aes = FALSE,
                size = text_size * 1.1,
                fontface = "bold",
                vjust = 0
            )
    }

    p <- p +
        ggdag::theme_dag() +
        # Add margins back (ggdag removes them)
        ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10, "mm")) +
        # prevent clipping of large nodes
        ggplot2::coord_cartesian(clip = "off") +
        # Expand axes to give space for nodes
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.2)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.2))

    if (length(unique(combined_dag_data$model_label)) > 1) {
        p <- p + ggplot2::facet_wrap(~model_label)
    }

    # --- Node Background Layer (drawn BEFORE edges so arrowheads show on top) ---
    # We draw filled squares/circles here so the arrow shaft is hidden inside the node,
    # but the arrowhead (at xend) is drawn ON TOP of this fill by the edge layers.
    node_data <- combined_dag_data |>
        dplyr::filter(!is.na(x) & !is.na(y)) |>
        dplyr::distinct(name, .keep_all = TRUE)

    p <- p +
        ggplot2::geom_point(
            data = node_data,
            ggplot2::aes(x = x, y = y, shape = type, fill = type),
            size   = uniform_node_size,
            color  = node_color,
            stroke = node_stroke,
            inherit.aes = FALSE
        ) +
        ggplot2::scale_shape_manual(
            values = c(Observed = 22, Latent = 21, Interaction = 23),
            guide = "none"
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                Observed    = node_fill,
                Latent      = node_fill,
                Interaction = "grey85"
            ),
            guide = "none"
        )

    if ("edge_type" %in% names(combined_dag_data)) {
        # Normalize bidirected edges to curve AWAY from the graph center.
        # Heuristic:
        # 1. Calculate Graph Centroid.
        # 2. For each edge, check if the "Right" side (curvature > 0) points towards or away from centroid.
        # 3. Swap direction if it points towards center.

        is_bidirected <- combined_dag_data$edge_type == "<->" &
            !is.na(combined_dag_data$edge_type)
        if (any(is_bidirected)) {
            # Calculate Centroid of the graph layout
            # Use unique node positions to avoid weighting by edge count
            unique_nodes <- unique(combined_dag_data[, c("name", "x", "y")])
            centroid_x <- mean(unique_nodes$x, na.rm = TRUE)
            centroid_y <- mean(unique_nodes$y, na.rm = TRUE)

            # Vectorized check
            # Current Vector P1 -> P2
            v_x <- combined_dag_data$xend - combined_dag_data$x
            v_y <- combined_dag_data$yend - combined_dag_data$y

            # Midpoint
            m_x <- (combined_dag_data$x + combined_dag_data$xend) / 2
            m_y <- (combined_dag_data$y + combined_dag_data$yend) / 2

            # Normal Vector pointing "Right" (Curvature direction)
            # R = (dy, -dx)
            r_x <- v_y
            r_y <- -v_x

            # Test Point: Midpoint + Normal
            t_x <- m_x + r_x
            t_y <- m_y + r_y

            # Distances to Centroid
            # dist_current: distance from (Midpoint + Right) to Centroid
            dist_current <- (t_x - centroid_x)^2 + (t_y - centroid_y)^2

            # dist_swapped: distance from (Midpoint - Right) to Centroid
            # (Midpoint - Right) coincides with the curve apex if we swapped direction
            t_swapped_x <- m_x - r_x
            t_swapped_y <- m_y - r_y
            dist_swapped <- (t_swapped_x - centroid_x)^2 +
                (t_swapped_y - centroid_y)^2

            # If swapped distance is GREATER (farther from center), we should swap
            # to make the curve point that way.
            needs_swap <- is_bidirected & (dist_swapped > dist_current)

            if (any(needs_swap, na.rm = TRUE)) {
                swap_idx <- which(needs_swap)

                # Temporary storage to enable swapping
                tmp_name <- combined_dag_data$name[swap_idx]
                tmp_x <- combined_dag_data$x[swap_idx]
                tmp_y <- combined_dag_data$y[swap_idx]

                combined_dag_data$name[swap_idx] <- combined_dag_data$to[
                    swap_idx
                ]
                combined_dag_data$x[swap_idx] <- combined_dag_data$xend[
                    swap_idx
                ]
                combined_dag_data$y[swap_idx] <- combined_dag_data$yend[
                    swap_idx
                ]

                combined_dag_data$to[swap_idx] <- tmp_name
                combined_dag_data$xend[swap_idx] <- tmp_x
                combined_dag_data$yend[swap_idx] <- tmp_y
            }
        }

        # Handle Bidirected Edges (Deduplicate A<->B vs B<->A)
        bidirected_edges <- combined_dag_data |>
            dplyr::filter(edge_type == "<->") |>
            dplyr::mutate(
                edge_id = paste(
                    pmin(name, to),
                    pmax(name, to),
                    sep = "_"
                )
            ) |>
            dplyr::distinct(edge_id, .keep_all = TRUE) |>
            dplyr::select(-edge_id)
        
        # Handle Directed Edges (Preserve All for Multi-Category Bundles)
        directed_edges <- combined_dag_data |>
            dplyr::filter(edge_type == "->")

        combined_dag_data <- dplyr::bind_rows(directed_edges, bidirected_edges)

        # 1. Directed Edges - using geom_dag_edges_arc which supports start/end caps.
        # By slicing data per curvature BEFORE building each layer, we avoid the
        # R lazy-eval bug where all layers would use the last curvature value.
        dir_edges <- directed_edges |>
            dplyr::mutate(
                dx    = xend - x,
                dy    = yend - y,
                mid_x = (x + xend) / 2 - curvature * dy / 2,
                mid_y = (y + yend) / 2 + curvature * dx / 2
            )

        if (nrow(dir_edges) > 0) {
            unique_curvatures <- sort(unique(dir_edges$curvature))
            edge_layers <- lapply(seq_along(unique_curvatures), function(i) {
                cv <- unique_curvatures[i]
                ld <- dir_edges[abs(dir_edges$curvature - cv) < 1e-6, , drop = FALSE]
                force(cv); force(ld)
                ggdag::geom_dag_edges_arc(
                    data       = ld,
                    mapping    = ggplot2::aes(
                        edge_width  = weight_abs,
                        edge_colour = significant,
                        label       = edge_label
                    ),
                    curvature   = cv,
                    angle_calc  = "along",
                    label_dodge = ggplot2::unit(3, "mm"),
                    start_cap   = cap_size,
                    end_cap     = cap_size
                )
            })
            p <- p + edge_layers
        }

        # 2. Bidirected Edges (grey, double arrow)
        if (nrow(bidirected_edges) > 0) {
            p <- p +
                ggdag::geom_dag_edges_arc(
                    data        = bidirected_edges,
                    mapping     = ggplot2::aes(label = edge_label),
                    curvature   = 0.5,
                    edge_width  = 0.3,
                    edge_colour = "grey60",
                    angle_calc  = "along",
                    label_dodge = ggplot2::unit(3, "mm"),
                    arrow       = ggplot2::arrow(
                        length = ggplot2::unit(2.5, "mm"),
                        type   = "closed",
                        ends   = "both"
                    ),
                    start_cap   = cap_size,
                    end_cap     = cap_size
                )
        }

        # Scales (ggraph edge aesthetics)
        p <- p +
            ggraph::scale_edge_colour_manual(
                values = c(
                    "pos"     = "firebrick",
                    "neg"     = "dodgerblue",
                    "sig"     = "black",
                    "ns"      = "grey70",
                    "default" = "black"
                ),
                guide = "none"
            ) +
            ggraph::scale_edge_width_continuous(
                range = edge_width_range,
                guide = "none"
            )
    } else {
        p <- p +
            ggdag::geom_dag_edges(
                start_cap = cap_size,
                end_cap = cap_size
            )
    }

    # --- Node Text Labels (Final Top Layer) ---
    p <- p +
        ggdag::geom_dag_text(
            data = node_data,
            ggplot2::aes(x = x, y = y, label = label_display),
            size   = text_size,
            colour = "black",
            inherit.aes = FALSE
        )

    return(p)
}


#' Convert Equations List to DAGitty String
#'
#' Implements the Attia, Holliday & Oldmeadow (2022) IDAG convention:
#' interaction terms (e.g. \code{BM:M}) and \code{I()} terms are represented as
#' explicit intermediate nodes rather than collapsed to direct component->response
#' edges.
#'
#' @param equations List of formulas
#' @param induced_cors List of character vectors (pairs) for bidirected edges
#' @param family Optional named character vector of family distributions
#' @return A named list:
#'   \itemize{
#'     \item \code{dag_string} — dagitty-compatible DAG string
#'     \item \code{interaction_nodes} — named list: internal_name -> display label
#'   }
#' @references
#'   Attia, J., Holliday, E., & Oldmeadow, C. (2022). A proposal for capturing
#'   interaction and effect modification using DAGs.
#'   \emph{International Journal of Epidemiology}, 51(4), 1047--1053.
#' @noRd
equations_to_dag_string <- function(
    equations,
    induced_cors = NULL,
    family = NULL,
    poly_terms = NULL, # list from get_all_polynomial_terms; used for fitted models
    # where equations are already expanded (I(age^2) -> age_pow2)
    collapse_expanded = FALSE
) {
    edges <- c()
    interaction_nodes <- list() # internal_name -> display label (e.g. "BM\u00d7M")

    # Build lookup: internal_name -> poly term info, for reconstructing diamond
    # nodes when equations are already expanded (fitted model path).
    poly_lookup <- list()
    if (!is.null(poly_terms)) {
        for (pt in poly_terms) {
            poly_lookup[[pt$internal_name]] <- pt
        }
    }

    occ_vars <- c()
    if (!is.null(family)) {
        occ_vars <- names(family)[family == "occupancy"]
    }

    # Helper: make a dagitty-safe node name from a term string
    # make_internal_name must match sanitize_term_name() from deterministic_nodes.R
    # so that node names in the plot align with JAGS parameter names like beta_weight_g_age_pow2.
    make_internal_name <- function(term) sanitize_term_name(term)

    # Helper: human-readable display label
    # For interactions: BM:M   -> BM×M
    # For I() powers:   I(age^2) -> age²  I(x^3) -> x³
    # For other I():    I(x+y)  -> I(x+y)  (keep as-is)
    make_display_label <- function(term) {
        if (grepl("^I\\(", term)) {
            # Check for simple power pattern: I(var^N)
            m <- regmatches(
                term,
                regexpr("^I\\(([a-zA-Z_][a-zA-Z0-9_]*)\\^([0-9]+)\\)$", term)
            )
            if (length(m) > 0) {
                inner <- sub("^I\\((.*)\\)$", "\\1", term)
                base <- sub("\\^.*$", "", inner)
                exp_n <- sub("^.*\\^", "", inner)
                superscripts <- c(
                    "\u2070",
                    "\u00b9",
                    "\u00b2",
                    "\u00b3",
                    "\u2074",
                    "\u2075",
                    "\u2076",
                    "\u2077",
                    "\u2078",
                    "\u2079"
                )
                n <- as.integer(exp_n)
                if (!is.na(n) && n >= 0 && n <= 9) {
                    return(paste0(base, superscripts[n + 1]))
                }
            }
            return(term) # fallback: keep as-is for complex I() expressions
        }
        gsub(":", "\u00d7", term) # "BM:M" -> "BM\u00d7M"
    }

    for (eq in equations) {
        resp <- all.vars(eq)[1]
        trm_lbls <- attr(terms(eq), "term.labels")

        # Optional: Collapse expanded names back to base names
        if (collapse_expanded) {
           resp <- gsub("(_[A-Za-z0-9]+|\\[[0-9]+\\])$", "", resp)
        }

        actual_resp <- if (resp %in% occ_vars) paste0("psi_", resp) else resp

        if (length(trm_lbls) == 0) {
            next
        } # intercept-only, nothing to draw

        # Detect pure deterministic declaration: entire RHS is a single I() call.
        # e.g.  AgeClass ~ I(0 * (age < 0.02) + 1 * (age >= 0.02))
        # Mark the LHS itself as the deterministic/interaction node and draw edges
        # from the component variables directly — avoids a giant intermediate node.
        is_pure_det <- length(trm_lbls) == 1 && grepl("^I\\(", trm_lbls[1])

        if (is_pure_det) {
            term <- trm_lbls[1]
            interaction_nodes[[actual_resp]] <- actual_resp # LHS is the det. node
            components <- all.vars(stats::as.formula(paste("~", term)))
            for (comp in components) {
                actual_comp <- if (comp %in% occ_vars) {
                    paste0("psi_", comp)
                } else {
                    comp
                }
                edges <- c(edges, paste(actual_resp, "<-", actual_comp))
            }
            next
        }

        for (term in trm_lbls) {
            # Skip random effects terms (e.g. 1 | year) which dagitty can't parse
            if (grepl("|", term, fixed = TRUE)) next
            
            # Optional: Collapse expanded names back to base names
            clean_term <- term
            if (collapse_expanded) {
               clean_term <- gsub("(_[A-Za-z0-9]+|\\[[0-9]+\\])$", "", term)
            }
            
            # Avoid self-loops if resp and term matched the same base name
            if (clean_term == resp) next

            is_interaction <- grepl(":", clean_term, fixed = TRUE) &&
                !grepl("^I\\(", clean_term)
            is_I_call <- grepl("^I\\(", clean_term)

            if (is_interaction || is_I_call) {
                # --- Deterministic / interaction node (within a mixed equation) ---
                # X:Y interactions and I() polynomial terms get a diamond node.
                iname <- make_internal_name(clean_term)
                d_label <- make_display_label(clean_term)
                interaction_nodes[[iname]] <- d_label

                # Component variables: split on : for interactions, all.vars for I()
                if (is_interaction) {
                    components <- strsplit(clean_term, ":", fixed = TRUE)[[1]]
                } else {
                    components <- all.vars(stats::as.formula(paste("~", clean_term)))
                }

                for (comp in components) {
                    actual_comp <- if (comp %in% occ_vars) {
                        paste0("psi_", comp)
                    } else {
                        comp
                    }
                    edges <- c(edges, paste(iname, "<-", actual_comp))
                }
                # Route through the deterministic node: response <- iname
                # For X:Y interactions, components may not appear as separate
                # explicit predictors, so also add comp -> resp edges.
                # For I() polynomial terms (e.g. I(age^2)), the base variable
                # (age) is typically an explicit separate term in the same equation,
                # so comp -> resp is already added by the regular-predictor branch.
                # Adding it again here creates a redundant hidden edge that
                # collides with the direct age -> weight_g arrow in the layout.
                edges <- c(edges, paste(actual_resp, "<-", iname))
                if (is_interaction) {
                    for (comp in components) {
                        actual_comp <- if (comp %in% occ_vars) {
                            paste0("psi_", comp)
                        } else {
                            comp
                        }
                        edges <- c(edges, paste(actual_resp, "<-", actual_comp))
                    }
                }
            } else {
                # --- Regular predictor ---
                # If we are collapsing, we skip the expansion-specific diamond nodes
                # and just draw an edge from the base variable to the response.
                if (collapse_expanded && term %in% names(poly_lookup)) {
                   # Skip the diamond node, handled by 'clean_term' falling through to Standard path
                } else if (!collapse_expanded && term %in% names(poly_lookup)) {
                    pt <- poly_lookup[[term]]
                    iname <- pt$internal_name # e.g. age_pow2
                    # ... (rest of diamond logic) ...
                    superscripts <- c("\u2070","\u00b9","\u00b2","\u00b3","\u2074","\u2075","\u2076","\u2077","\u2078","\u2079")
                    n <- as.integer(pt$power)
                    d_label <- if (!is.na(n) && n >= 0 && n <= 9) paste0(pt$base_var, superscripts[n+1]) else iname
                    interaction_nodes[[iname]] <- d_label
                    base_actual <- if (pt$base_var %in% occ_vars) paste0("psi_", pt$base_var) else pt$base_var
                    edges <- c(edges, paste(iname, "<-", base_actual))
                    edges <- c(edges, paste(actual_resp, "<-", iname))
                } else {
                    actual_pred <- if (clean_term %in% occ_vars) {
                        paste0("psi_", clean_term)
                    } else {
                        clean_term
                    }
                    edges <- c(edges, paste(actual_resp, "<-", actual_pred))
                }
            }
        }
    }

    # Occupancy p_ nodes (visual grouping only)
    for (v in occ_vars) {
        # p_node existence is implied by its edges elsewhere
    }

    if (!is.null(induced_cors)) {
        for (pair in induced_cors) {
            if (length(pair) == 2) {
                p1 <- if (pair[1] %in% occ_vars) {
                    paste0("psi_", pair[1])
                } else {
                    pair[1]
                }
                p2 <- if (pair[2] %in% occ_vars) {
                    paste0("psi_", pair[2])
                } else {
                    pair[2]
                }
                edges <- c(edges, paste(p1, "<->", p2))
            }
        }
    }

    list(
        dag_string = paste("dag {", paste(unique(edges), collapse = "; "), "}"),
        interaction_nodes = interaction_nodes
    )
}
