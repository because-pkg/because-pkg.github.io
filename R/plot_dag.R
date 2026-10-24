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
#'   Override node positions entirely with the \code{coords} argument. Use \code{"multiscale"}
#'   to automatically arrange nodes in vertical tiers based on their hierarchical level.
#' @param latent Character vector of latent variable names. Overrides the model's
#'   latent variables if provided. If \code{NULL} (default) and plotting from formulas,
#'   latent variables are **automatically detected** if they follow the SEM
#'   naming convention (e.g., \code{L1}, \code{L2}, \code{Latent1}, \code{lat_foo})
#'   and only appear on the RHS of equations.
#' @param node_size Size of the nodes (default 14).
#' @param node_color Color of the node border (default "black").
#' @param node_fill Color of the node interior (default "white").
#' @param node_stroke Thickness of the node border (default 1.5).
#' @param text_size Size of the node labels (default 3.5).
#' @param edge_label_size Size of the edge coefficient labels (default 3). Reduce for crowded plots (e.g., 2 or 2.5).
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
#' @param type Character; either "raw" (default) or "marginal" effects.
#' @param multinomial_probabilities Logical; if TRUE, expands multinomial categories.
#'
#' @return A `ggplot` object that can be further customized with standard ggplot2 functions (e.g., `+ theme_...()`, `+ ggtitle(...)`).
#'
#' @details
#' \strong{Interaction and Deterministic Nodes:}
#' Interaction terms (e.g. \code{BM:M}) and \code{I()} transformations are
#' rendered as **explicit intermediate nodes** (grey diamonds), following the
#' Interaction DAG (IDAG) convention of Attia, Holliday \& Oldmeadow (2022).
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
#' Attia, J., Holliday, E., \& Oldmeadow, C. (2022). A proposal for capturing
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
#' @export
plot_dag <- function(
    x,
    layout = "kk",
    latent = NULL,
    node_size = 12,
    node_color = "black",
    node_fill = "white",
    node_stroke = 1.2,
    text_size = 3.5,
    edge_label_size = 3,
    edge_width_range = c(0.5, 2),
    edge_color_scheme = c("directional", "binary", "monochrome"),
    show_coefficients = TRUE,
    coords = NULL,
    family = NULL,
    type = c("raw", "marginal"),
    multinomial_probabilities = TRUE
) {
    edge_color_scheme <- match.arg(edge_color_scheme)
    type              <- match.arg(type)

    # Check dependencies
    if (
        !requireNamespace("dagitty", quietly = TRUE) ||
            !requireNamespace("ggdag", quietly = TRUE) ||
            !requireNamespace("ggraph", quietly = TRUE) ||
            !requireNamespace("ggplot2", quietly = TRUE) ||
            !requireNamespace("dplyr", quietly = TRUE)
    ) {
        stop("Packages 'dagitty', 'ggdag', 'ggraph', 'ggplot2', and 'dplyr' are required for plot_dag.")
    }

    # Normalize input to a named list of model objects / formula lists
    if (inherits(x, "because") || inherits(x, "list") && all(sapply(x, inherits, "formula"))) {
        x <- list(Model = x)
    }

    # --- Build tidy data frame for each model ---
    combined_dag_data <- NULL

    for (i in seq_along(x)) {
        obj   <- x[[i]]
        label <- names(x)[i]
        if (is.null(label)) label <- paste("Model", i)

        dag_data <- build_dag_data(
            obj, label,
            latent                   = latent,
            family                   = family,
            layout                   = layout,
            coords                   = coords,
            node_size                = node_size,
            text_size                = text_size,
            edge_label_size          = edge_label_size,
            edge_color_scheme        = edge_color_scheme,
            type                     = type,
            multinomial_probabilities = multinomial_probabilities,
            show_coefficients        = show_coefficients
        )

        if (is.null(combined_dag_data)) {
            combined_dag_data <- dag_data
        } else {
            combined_dag_data <- rbind(combined_dag_data, dag_data)
        }
    }

    # --- Compound bounding boxes (imperfect detection / state-space groups) ---
    compound_boxes <- NULL
    if (
        !is.null(combined_dag_data) &&
            "compound_group" %in% names(combined_dag_data) &&
            any(!is.na(combined_dag_data$compound_group))
    ) {
        if (any(!is.na(combined_dag_data$x))) {
            compound_boxes <- combined_dag_data |>
                dplyr::filter(!is.na(compound_group)) |>
                dplyr::group_by(model_label, compound_group) |>
                dplyr::summarize(
                    xmin = min(x) - 0.25,
                    xmax = max(x) + 0.25,
                    ymin = min(y) - 0.25,
                    ymax = max(y) + 0.25,
                    .groups = "drop"
                )
        }
    }

    # --- Finalize node type factor ---
    combined_dag_data$type <- factor(
        combined_dag_data$type,
        levels = c("Observed", "Latent", "Interaction")
    )

    # --- Uniform node size ---
    if (is.null(combined_dag_data$final_node_size)) {
        uniform_node_size <- node_size
    } else {
        uniform_node_size <- max(combined_dag_data$final_node_size, na.rm = TRUE)
        if (!is.finite(uniform_node_size)) uniform_node_size <- node_size
    }

    cap_size <- ggraph::circle(uniform_node_size * 1.15, "pt")

    # --- ggplot assembly ---
    p <- ggplot2::ggplot(
        combined_dag_data,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
    )

    # 0. Compound group bounding boxes (background layer)
    if (!is.null(compound_boxes)) {
        p <- p +
            ggplot2::geom_rect(
                data = compound_boxes,
                ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE,
                fill = "grey93", color = "grey80",
                alpha = 0.5, linewidth = 0.4, lty = "dashed"
            ) +
            ggplot2::geom_text(
                data = compound_boxes,
                ggplot2::aes(x = (xmin + xmax) / 2, y = ymax + 0.05, label = compound_group),
                inherit.aes = FALSE,
                size = text_size * 1.1, fontface = "bold", vjust = 0
            )
    }

    p <- p +
        ggdag::theme_dag() +
        ggplot2::theme(plot.margin = ggplot2::margin(10, 10, 10, 10, "mm")) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.25)) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.25))

    if (length(unique(combined_dag_data$model_label)) > 1) {
        p <- p + ggplot2::facet_wrap(~model_label)
    }

    # Node background layer (drawn BEFORE edges so arrowheads show on top)
    node_data <- combined_dag_data |>
        dplyr::filter(!is.na(x) & !is.na(y)) |>
        dplyr::distinct(name, .keep_all = TRUE)

    p <- p +
        ggplot2::geom_point(
            data   = node_data,
            ggplot2::aes(x = x, y = y, shape = type, fill = type),
            size   = uniform_node_size, color = node_color,
            stroke = node_stroke, inherit.aes = FALSE
        ) +
        ggplot2::scale_shape_manual(
            values = c(Observed = 22, Latent = 21, Interaction = 23),
            guide  = "none"
        ) +
        ggplot2::scale_fill_manual(
            values = c(Observed = node_fill, Latent = node_fill, Interaction = "grey85"),
            guide  = "none"
        )

    # Edge layers
    if ("edge_type" %in% names(combined_dag_data)) {
        # Orient bidirected edges away from graph centroid
        is_bidirected <- combined_dag_data$edge_type == "<->" & !is.na(combined_dag_data$edge_type)
        if (any(is_bidirected)) {
            unique_nodes <- unique(combined_dag_data[, c("name", "x", "y")])
            centroid_x   <- mean(unique_nodes$x, na.rm = TRUE)
            centroid_y   <- mean(unique_nodes$y, na.rm = TRUE)

            v_x <- combined_dag_data$xend - combined_dag_data$x
            v_y <- combined_dag_data$yend - combined_dag_data$y
            m_x <- (combined_dag_data$x + combined_dag_data$xend) / 2
            m_y <- (combined_dag_data$y + combined_dag_data$yend) / 2
            r_x <- v_y;  r_y <- -v_x
            t_x <- m_x + r_x;  t_y <- m_y + r_y

            dist_current <- (t_x - centroid_x)^2 + (t_y - centroid_y)^2
            t_swapped_x  <- m_x - r_x;  t_swapped_y <- m_y - r_y
            dist_swapped <- (t_swapped_x - centroid_x)^2 + (t_swapped_y - centroid_y)^2
            needs_swap   <- is_bidirected & (dist_swapped > dist_current)

            if (any(needs_swap, na.rm = TRUE)) {
                swap_idx <- which(needs_swap)
                tmp_name <- combined_dag_data$name[swap_idx]
                tmp_x    <- combined_dag_data$x[swap_idx]
                tmp_y    <- combined_dag_data$y[swap_idx]
                combined_dag_data$name[swap_idx]  <- combined_dag_data$to[swap_idx]
                combined_dag_data$x[swap_idx]     <- combined_dag_data$xend[swap_idx]
                combined_dag_data$y[swap_idx]     <- combined_dag_data$yend[swap_idx]
                combined_dag_data$to[swap_idx]    <- tmp_name
                combined_dag_data$xend[swap_idx]  <- tmp_x
                combined_dag_data$yend[swap_idx]  <- tmp_y
            }
        }

        # Deduplicate bidirected edges (A<->B vs B<->A)
        bidirected_edges <- combined_dag_data |>
            dplyr::filter(edge_type == "<->") |>
            dplyr::mutate(edge_id = paste(pmin(name, to), pmax(name, to), sep = "_")) |>
            dplyr::distinct(edge_id, .keep_all = TRUE) |>
            dplyr::select(-edge_id)

        directed_edges   <- combined_dag_data |>
            dplyr::filter(edge_type == "->")

        combined_dag_data <- dplyr::bind_rows(directed_edges, bidirected_edges)

        # Directed edges — one layer per unique curvature (avoids R lazy-eval scoping bug)
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
                    data      = ld,
                    mapping   = ggplot2::aes(
                        edge_width  = weight_abs,
                        edge_colour = significant,
                        label       = edge_label
                    ),
                    curvature   = cv,
                    angle_calc  = "along",
                    label_dodge = ggplot2::unit(3, "mm"),
                    label_size  = edge_label_size,
                    start_cap   = cap_size,
                    end_cap     = cap_size
                )
            })
            p <- p + edge_layers
        }

        # Bidirected edges (grey, double arrow)
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
                    label_size  = edge_label_size,
                    arrow       = ggplot2::arrow(
                        length = ggplot2::unit(2.5, "mm"),
                        type   = "closed", ends = "both"
                    ),
                    start_cap   = cap_size,
                    end_cap     = cap_size
                )
        }

        p <- p +
            ggraph::scale_edge_colour_manual(
                values = c(pos = "dodgerblue", neg = "firebrick", sig = "black", ns = "grey70", default = "black"),
                guide  = "none"
            ) +
            ggraph::scale_edge_width_continuous(range = edge_width_range, guide = "none")

    } else {
        p <- p + ggdag::geom_dag_edges(start_cap = cap_size, end_cap = cap_size)
    }

    # Node text labels (top layer)
    p <- p +
        ggdag::geom_dag_text(
            data       = node_data,
            ggplot2::aes(x = x, y = y, label = label_display),
            size       = text_size,
            colour     = "black",
            inherit.aes = FALSE
        )

    return(p)
}
