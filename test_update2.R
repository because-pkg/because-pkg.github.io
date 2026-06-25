devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")
eqs <- list(Y ~ X)
dat <- data.frame(Y = rnorm(100), X = rnorm(100))
fit <- because(eqs, data = dat, n.iter = 2000)

cat("\nUpdating model...\n")
fit_new <- try(update(fit, engine = "numpyro", n.iter=100))
if(inherits(fit_new, "try-error")) {
  reticulate::py_last_error()
}
