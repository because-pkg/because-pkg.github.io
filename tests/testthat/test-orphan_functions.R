library(testthat)
library(because)

# Create a mock because object for testing functions that require it
create_mock_because <- function() {
    skip_if_not_installed("rjags")
    set.seed(123)
    df <- data.frame(X = rnorm(20), Y = rnorm(20))
    # Quick jags run without NA to allow pp_check to work
    fit <- because(list(Y ~ X), data = df, engine = "jags", n.iter = 100, n.burnin = 10, quiet = TRUE, monitor = "all")
    return(fit)
}

test_that("Core modeling orphan functions execute without errors", {
    df <- data.frame(X = 1:10, Y = 1:10)
    
    # because_format_data (general id_col + structure and backward-compatible species_col)
    df_format <- data.frame(
        ID = rep(c("A", "B"), each = 5),
        Const = rep(c(10, 20), each = 5),
        Y = 1:10
    )
    fmt <- because_format_data(df_format, id_col = "ID")
    expect_type(fmt, "list")
    expect_true(is.vector(fmt$Const))
    expect_true(is.matrix(fmt$Y))
    expect_equal(dim(fmt$Y), c(2, 5))

    # Aligning with a matrix structure
    struct_mat <- diag(2)
    rownames(struct_mat) <- colnames(struct_mat) <- c("B", "A")
    fmt_struct <- because_format_data(df_format, id_col = "ID", structure = struct_mat)
    expect_equal(rownames(fmt_struct$Y), c("B", "A"))
    expect_equal(names(fmt_struct$Const), c("B", "A"))

    # Backward compatibility with species_col = "SP"
    df_legacy <- data.frame(SP = rep(c("A", "B"), each = 3), Y = 1:6)
    fmt_legacy <- because_format_data(df_legacy, species_col = "SP")
    expect_type(fmt_legacy, "list")
    
    # because_structure (with and without name)
    custom_struct <- because_structure(
        precision_fn = function(n) diag(n)
    )
    expect_type(custom_struct, "closure")
    custom_struct_named <- because_structure(
        name = "test_struct",
        precision_fn = function(n) diag(n)
    )
    expect_type(custom_struct_named, "closure")
    
    # because_family (with and without name)
    custom_fam <- because_family(
        jags_likelihood = "{response}[{i}] ~ dnorm({mu}[{i}], {tau})"
    )
    expect_type(custom_fam, "closure")
    custom_fam_named <- because_family(
        name = "test_fam",
        jags_likelihood = "{response}[{i}] ~ dnorm({mu}[{i}], {tau})"
    )
    expect_type(custom_fam_named, "closure")
})

test_that("Post-estimation orphan functions execute without errors", {
    fit <- create_mock_because()
    
    # marginal_effects
    me <- marginal_effects(fit)
    expect_s3_class(me, "data.frame")
    
    # because_mediation
    expect_error(because_mediation(fit, "X", "Y"), NA)
    
    # posterior_predict
    expect_error(posterior_predict(fit), NA)
})

test_that("Plotting and diagnostics orphan functions execute without errors", {
    skip_if_not_installed("rjags")
    skip_if_not_installed("ggplot2")
    fit <- create_mock_because()
    expect_s3_class(plot_coef(fit), "ggplot")
    expect_error(plot_posterior(fit, "beta_Y_X"), NA)
    expect_s3_class(pp_check(fit, resp = "Y"), "ggplot")
})
