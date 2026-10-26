test_that("because_dsep identifies correct basis set for simple DAG", {
    # A -> B -> C
    equations <- list(
        B ~ A,
        C ~ B
    )

    result <- because_dsep(equations)

    # Should imply A _||_ C | B
    expect_type(result, "list")
    expect_true(length(result) > 0)

    # Check structure of result (it returns a list of formulas)
    expect_s3_class(result[[1]], "formula")
})

test_that("because_dsep returns empty list for saturated model", {
    # A -> B -> C, A -> C
    equations <- list(
        B ~ A,
        C ~ A + B
    )

    result <- because_dsep(equations)

    expect_length(result, 0)
})

test_that("because_dsep handles latent variables", {
    # L -> X, L -> Y (L is latent)
    # Should imply X _||_ Y is FALSE (they are correlated)
    # But d-sep on MAG should show they are connected by bidirected edge

    equations <- list(
        X ~ L,
        Y ~ L
    )

    # This should return the induced correlations
    result <- because_dsep(equations, latent = "L")

    expect_type(result, "list")
    # Check that result contains at least tests and correlations
    expect_true(all(c("tests", "correlations") %in% names(result)))

    # No d-sep tests (saturated in MAG sense for these 2 vars?)
    # Actually X <-> Y means no independence

    # Check correlations
    expect_true(length(result$correlations) > 0)
    expect_equal(sort(result$correlations[[1]]), c("X", "Y"))
})
