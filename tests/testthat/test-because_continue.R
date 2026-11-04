test_that("because_continue() works with sequential JAGS", {
  skip_if_not_installed("rjags")
  set.seed(42)
  df <- data.frame(X = rnorm(30), Y = rnorm(30))

  fit <- because(
    equations = list(Y ~ X),
    data      = df,
    engine    = "jags",
    n.chains  = 2,
    n.iter    = 50,
    n.burnin  = 0,
    n.thin    = 1,
    parallel  = FALSE,
    quiet     = TRUE
  )

  expect_s3_class(fit, "because")
  n_orig <- nrow(as.matrix(fit$samples[[1]]))
  expect_equal(n_orig, 50)

  # Continue sampling
  fit_cont <- because_continue(fit, n.iter = 50, quiet = TRUE)
  expect_s3_class(fit_cont, "because")
  expect_equal(nrow(as.matrix(fit_cont$samples[[1]])), 100)
  expect_true("beta_Y_X" %in% rownames(fit_cont$summary$statistics))
})

test_that("because_continue() works with parallel JAGS rebuild when model object is saved/restored", {
  skip_if_not_installed("rjags")
  set.seed(42)
  df <- data.frame(X = rnorm(30), Y = rnorm(30))

  fit <- because(
    equations = list(Y ~ X),
    data      = df,
    engine    = "jags",
    n.chains  = 2,
    n.iter    = 50,
    n.burnin  = 0,
    n.thin    = 1,
    parallel  = TRUE,
    n.cores   = 2,
    quiet     = TRUE
  )

  expect_s3_class(fit, "because")

  # Simulate saving and reloading RDS (which invalidates live external pointers)
  tmp_rds <- tempfile(fileext = ".rds")
  saveRDS(fit, tmp_rds)
  on.exit(unlink(tmp_rds), add = TRUE)
  reloaded <- readRDS(tmp_rds)

  # Rebuild in parallel
  fit_cont <- because_continue(reloaded, n.iter = 50, parallel = TRUE, n.cores = 2, quiet = TRUE)
  expect_s3_class(fit_cont, "because")
  expect_equal(nrow(as.matrix(fit_cont$samples[[1]])), 100)
  expect_true("beta_Y_X" %in% rownames(fit_cont$summary$statistics))
})
