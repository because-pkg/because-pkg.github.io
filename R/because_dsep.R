#' Extract d-separation statements from a structural equation model
#'
#' This function takes a set of structural equations defining a causal model
#' and returns the conditional independence statements (d-separation or m-separation tests)
#' implied by the model structure. If latent variables are specified, the function
#' uses the MAG (Mixed Acyclic Graph) approach by Shipley and Douma (2021)
#' to account for unmeasured latent variables.
#'
#' @param equations A list of model formulas (one per structural equation),
#'   e.g., \code{list(Y ~ X1 + X2, Z ~ Y)}.
#' @param latent Optional character vector of latent (unmeasured) variable names.
#'   If provided, the function converts the DAG to a MAG and returns m-separation tests.
#'
#' @return If \code{latent} is NULL, returns a list of formulas representing
#'   conditional independence tests. If \code{latent} is specified, returns a list with:
#'   \itemize{
#'     \item \code{tests}: List of m-separation test formulas
#'     \item \code{correlations}: List of variable pairs with induced correlations
#'   }
#'
#' @details
#' The function implements the basis set approach to d-separation testing
#' (Shipley 2000, 2009, 2016). For standard DAGs without latent variables, it identifies
#' pairs of non-adjacent variables and creates conditional independence tests.
#'
#' When latent variables are specified, the function uses the DAG-to-MAG conversion
#' (Shipley & Douma 2021) to identify m-separation statements and induced correlations
#' among observed variables that arise from shared latent common causes.
#'
#' Deterministic nodes (interaction terms such as \code{A:B}, and arithmetic
#' transformations such as \code{I(A^2)}) are kept as **explicit intermediate
#' nodes** in the DAG, following the D-separation extension of
#' Geiger, Verma & Pearl (1990).  This ensures that the basis set includes
#' independence tests that condition on the deterministic term itself
#' (e.g. \eqn{TL \perp BM \mid \{BM{:}M\}}), which would be silently dropped
#' if the interaction were collapsed to its component variables.
#'
#' @references
#' Geiger, D., Verma, T., & Pearl, J. (1990). Identifying independence in
#' Bayesian Networks. \emph{Networks}, 20(5), 507--534.
#'
#' Shipley, B. (2000). A new inferential test for path models based on
#' directed acyclic graphs. Structural Equation Modeling, 7(2), 206-218.
#'
#' Shipley, B. (2009). Confirmatory path analysis in a generalized multilevel
#' context. Ecology, 90(2), 363-368.
#'
#' Shipley, B. (2016). Cause and Correlation in Biology (2nd ed.).
#' Cambridge University Press.
#'
#' Shipley, B., & Douma, J. C. (2021). Testing Piecewise Structural Equations
#' Models in the Presence of Latent Variables and Including Correlated Errors.
#' Structural Equation Modeling: A Multidisciplinary Journal, 28(4), 582-589.
#' https://doi.org/10.1080/10705511.2020.1871355
#'

#'
#' @examples
#' # Standard DAG
#' equations <- list(LS ~ BM, NL ~ BM + RS, DD ~ NL)
#' ind_tests <- because_dsep(equations)
#'
#' # With latent variable
#' equations_latent <- list(X ~ Quality, Y ~ Quality)
#' result <- because_dsep(equations_latent, latent = "Quality")
#' # result$tests: m-separation tests
#' # result$correlations: induced correlation between X and Y
#' @param poly_terms Internal list of polynomial terms.
#' @param categorical_vars Character vector of categorical variable names.
#' @param family Named character vector of family/distribution for response variables.
#' @param quiet Logical; if FALSE (default), print the basis set and MAG structure.
#'   If TRUE, suppress informational output.
#' @param random_terms Optional list of random effects (group, type) parsed from equations.
#' @param hierarchical_info Internal argument used to pass data hierarchy information
#'   (levels, grouping variables) for multiscale d-separation
#'   tests (following Shipley 2009).
#' @export
#' @importFrom stats formula terms as.formula
because_dsep <- function(
  equations,
  latent = NULL,
  random_terms = list(),
  hierarchical_info = NULL,
  poly_terms = NULL,
  categorical_vars = NULL,
  family = NULL,
  quiet = FALSE
) {
  # Support passing a fitted 'because' object directly
  if (inherits(equations, "because")) {
    obj <- equations
    equations <- obj$equations
    if (is.null(latent) && !is.null(obj$latent)) latent <- obj$latent
    if (is.null(hierarchical_info) && !is.null(obj$hierarchical_info)) hierarchical_info <- obj$hierarchical_info
  }

  # If no latents, use standard DAG d-separation
  if (is.null(latent)) {
    return(dsep_standard(
      equations,
      random_terms = random_terms,
      hierarchical_info = hierarchical_info,
      poly_terms = poly_terms,
      categorical_vars = categorical_vars,
      family = family,
      quiet = quiet
    ))
  }

  # With latents: use MAG m-separation
  return(dsep_with_latents(
    equations,
    latent,
    random_terms = random_terms,
    hierarchical_info = hierarchical_info,
    poly_terms = poly_terms,
    categorical_vars = categorical_vars,
    family = family,
    quiet = quiet
  ))
}



#' plot_dsep
#'
#' Creates a caterpillar plot (point and whisker) of the regression coefficients
#' from all d-separation tests. A horizontal red line at zero helps visually
#' assess which independence claims are fulfilled (95% CI includes zero) or
#' violated (95% CI excludes zero).
#'
#' @param object A `because` object fitted with \code{dsep = TRUE}.
#' @param ... Additional arguments.
#' @param prob Numeric; probability mass for the credibility interval (default 0.95).
#'
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' # Plot results for a fitted model
#' plot_dsep(fit)
#' }
#' @rdname plot_dsep
#' @export
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_hline coord_flip labs theme_minimal theme
plot_dsep.because <- function(object, prob = 0.95, ...) {
  has_dsep <- !is.null(object$dsep) &&
    (isTRUE(object$dsep) ||
       (is.list(object$dsep) && !is.null(object$dsep$results)))
  if (!has_dsep) {
    stop("plot_dsep requires a 'because' object fitted with dsep = TRUE.")
  }

  # summary.because handles the complex renaming and dummy variable matching
  s <- summary(object, prob = prob)

  if (is.null(s$results) || nrow(s$results) == 0) {
    stop("No d-separation test results found in model object.")
  }

  res <- s$results

  # Ensure ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_dsep. Please install it.")
  }

  # For multinomial responses, multiple parameters share the same Test string
  # (one per category). Disambiguate by appending the specific Parameter name
  # so that each category gets its own row in the caterpillar plot.
  dup_tests <- duplicated(res$Test) | duplicated(res$Test, fromLast = TRUE)
  if (any(dup_tests)) {
    res$Test[dup_tests] <- paste0(res$Test[dup_tests], " (", res$Parameter[dup_tests], ")")
  }

  # Create Simplified labels for the plot: Response _||_ TestVar
  # This strips the conditioning set (e.g. | {X,Y,Z}) for better readability.
  # We match the specific formatting used in summary.because.
  res$Label <- gsub(" \\| \\{.*?\\}", "", res$Test)
  res$Label <- trimws(res$Label)

  # Symmetric axis: always show both sides of zero for fair visual comparison.
  # Use coord_flip(ylim=) rather than scale_y_continuous(limits=) to zoom
  # the coordinate system without clipping whiskers at the scale level.
  max_abs <- max(abs(c(res$LowerCI, res$UpperCI, res$Estimate)), na.rm = TRUE)
  axis_lim <- c(-max_abs, max_abs) * 1.05  # 5% padding

  # Create Plot
  p <- ggplot2::ggplot(
    res,
    ggplot2::aes(
      x = stats::reorder(Label, seq_len(nrow(res))),
      y = Estimate,
      ymin = LowerCI,
      ymax = UpperCI
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey50",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    ggplot2::geom_pointrange(linewidth = 0.8, size = 0.5, fatten = 3) +
    ggplot2::coord_flip(ylim = axis_lim) +
    ggplot2::labs(
      title = "d-separation Independence Tests",
      subtitle = paste0("Caterpillar plot of path coefficients with ", prob * 100, "% Bayesian Credibility Intervals"),
      x = "Conditional Independence Claim",
      y = "Estimated Beta (Effect Size)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 11),
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}

