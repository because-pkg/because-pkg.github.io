# Tests for Geiger et al. 1990 D-separation compliance with deterministic nodes
# Verifies that interaction terms (e.g. BM:M) are kept as explicit intermediate
# nodes in the DAG so that the full Geiger basis set is produced, including tests
# that condition on the interaction term itself.

library(testthat)

test_that("basis set for interaction DAG follows Geiger et al. 1990", {
    # DAG:  LS ~ BM, TL ~ BM:M, M ~ 1
    # After Geiger fix the explicit BM_x_M node in the DAG produces 6 tests:
    #   1. BM  _||_ M           |  {}
    #   2. BM  _||_ TL          | {BM:M}  <- NEW: Geiger-compliant test
    #   3. LS  _||_ M           | {BM}
    #   4. LS  _||_ BM:M        | {BM, M} (interaction node vs LS)
    #   5. LS  _||_ TL          | {BM, BM:M}
    #   6. M   _||_ TL          | {BM:M}  <- NEW: Geiger-compliant test
    # (exact order may vary by ggm version; we check key structural properties)

    equations <- list(LS ~ BM, TL ~ BM:M, M ~ 1)
    result <- because_dsep(equations, quiet = TRUE)

    expect_type(result, "list")

    # Pre-fix there were only 3 tests; Geiger compliance adds 3 more
    expect_gt(
        length(result),
        3,
        label = "Geiger-compliant basis set should have more than 3 tests"
    )

    # Helper: does a formula's RHS contain BM:M (or its internal name)?
    has_interaction_in_rhs <- function(f) {
        rhs_str <- as.character(f)[3]
        grepl("BM:M|BM_x_M", rhs_str)
    }

    # Geiger test 4: {BM, TL} pair conditioned on {BM:M}
    # (may appear as BM ~ TL + BM:M OR TL ~ BM + BM:M)
    geiger_test_bm <- any(vapply(
        result,
        function(f) {
            all_v <- all.vars(f)
            "BM" %in% all_v && "TL" %in% all_v && has_interaction_in_rhs(f)
        },
        logical(1)
    ))

    # Geiger test 5: {M, TL} pair conditioned on {BM:M}
    # (may appear as M ~ TL + BM:M OR TL ~ M + BM:M)
    geiger_test_m <- any(vapply(
        result,
        function(f) {
            resp <- as.character(f)[2]
            test_var <- attr(f, "test_var")
            pair <- sort(c(resp, test_var))
            identical(pair, c("M", "TL")) && has_interaction_in_rhs(f)
        },
        logical(1)
    ))

    expect_true(
        geiger_test_bm,
        label = "Basis set must include {BM,TL} test conditioned on BM:M (Geiger 1990 test 4)"
    )

    expect_true(
        geiger_test_m,
        label = "Basis set must include {M,TL} test conditioned on BM:M (Geiger 1990 test 5)"
    )
})

test_that("because_dsep still produces correct basis set for simple chain without interactions", {
    # A -> B -> C  =>  A _||_ C | {B}
    equations <- list(B ~ A, C ~ B)
    result <- because_dsep(equations, quiet = TRUE)

    expect_type(result, "list")
    expect_length(result, 1)

    # The single test should involve A, B, C
    all_v <- all.vars(result[[1]])
    expect_true(all(c("A", "B", "C") %in% all_v))
})

test_that("saturated DAG has empty basis set after Geiger fix", {
    # A -> B -> C, A -> C  (saturated: no missing edges)
    equations <- list(B ~ A, C ~ A + B)
    result <- because_dsep(equations, quiet = TRUE)
    expect_length(result, 0)
})
