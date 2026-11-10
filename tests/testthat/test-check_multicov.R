library(because)

testthat::test_that("check_multicov computes Rhat and Law of Total Variance decomposition correctly", {
  testthat::skip_if_not_installed("rjags")
  
  set.seed(42)
  N <- 12
  x <- rnorm(N)
  y <- 0.6 * x + rnorm(N, sd = 0.5)
  dat <- data.frame(ID = paste0("S", 1:N), x = x, y = y)

  coords1 <- cbind(runif(N), runif(N))
  coords2 <- cbind(runif(N) * 2, runif(N) * 2)
  m1 <- exp(-as.matrix(dist(coords1)))
  m2 <- exp(-as.matrix(dist(coords2)))
  rownames(m1) <- colnames(m1) <- dat$ID
  rownames(m2) <- colnames(m2) <- dat$ID

  # 1. Single collection
  fit_single <- because(
    equations = list(y ~ x),
    data = dat,
    structure = list(Kernel_Exp = m1, Kernel_Gaus = m2),
    id_col = "ID",
    engine = "jags",
    n.chains = 2,
    n.iter = 300,
    n.burnin = 100,
    quiet = TRUE
  )

  diag_single <- check_multicov(fit_single, plot = FALSE)
  testthat::expect_s3_class(diag_single, "check_multicov")
  testthat::expect_true(!is.null(diag_single$variance_decomposition))
  pct_w_name <- intersect(names(diag_single$variance_decomposition), c("Pct_within_tree", "Pct_within_struct", "Pct_MCMC"))[1]
  pct_s_name <- intersect(names(diag_single$variance_decomposition), c("Pct_between_trees", "Pct_between_structs", "Pct_Struct"))[1]
  testthat::expect_true(all(c("Parameter", "SD", pct_w_name, pct_s_name) %in% names(diag_single$variance_decomposition)))
  
  pct_sum_single <- diag_single$variance_decomposition[[pct_w_name]] + diag_single$variance_decomposition[[pct_s_name]]
  testthat::expect_true(all(abs(pct_sum_single - 100.0) <= 0.2))

  # 2. Nested multi-group collection
  coords3 <- cbind(runif(N) * 3, runif(N) * 3)
  coords4 <- cbind(runif(N) * 4, runif(N) * 4)
  m3 <- exp(-as.matrix(dist(coords3)))
  m4 <- exp(-as.matrix(dist(coords4)))
  rownames(m3) <- colnames(m3) <- dat$ID
  rownames(m4) <- colnames(m4) <- dat$ID

  fit_grouped <- because(
    equations = list(y ~ x),
    data = dat,
    structure = list(
      GroupA = list(Tree1 = m1, Tree2 = m2),
      GroupB = list(Tree3 = m3, Tree4 = m4)
    ),
    id_col = "ID",
    engine = "jags",
    n.chains = 2,
    n.iter = 300,
    n.burnin = 100,
    quiet = TRUE
  )

  diag_grouped <- check_multicov(fit_grouped, plot = FALSE)
  testthat::expect_s3_class(diag_grouped, "check_multicov")
  testthat::expect_true(!is.null(diag_grouped$variance_decomposition))
  pct_gw_name <- intersect(names(diag_grouped$variance_decomposition), c("Pct_within_tree", "Pct_within_struct", "Pct_MCMC"))[1]
  pct_gs_name <- intersect(names(diag_grouped$variance_decomposition), c("Pct_between_trees", "Pct_between_structs", "Pct_Struct"))[1]
  pct_gg_name <- intersect(names(diag_grouped$variance_decomposition), c("Pct_backbone", "Pct_group", "Pct_Group"))[1]
  testthat::expect_true(all(c("Parameter", "SD", pct_gw_name, pct_gs_name, pct_gg_name) %in% names(diag_grouped$variance_decomposition)))
  
  pct_sum_grouped <- diag_grouped$variance_decomposition[[pct_gw_name]] + 
                     diag_grouped$variance_decomposition[[pct_gs_name]] + 
                     diag_grouped$variance_decomposition[[pct_gg_name]]
  testthat::expect_true(all(abs(pct_sum_grouped - 100.0) <= 0.2))

  # Test plot method does not fail
  tf <- tempfile(fileext = ".pdf")
  pdf(tf)
  testthat::expect_silent(plot(diag_grouped))
  dev.off()
  if (file.exists(tf)) unlink(tf)
})
