test_that("get_family_default_rope returns expected principled bounds", {
    # Gaussian / continuous
    expect_equal(get_family_default_rope("gaussian"), c(-0.10, 0.10))
    expect_equal(get_family_default_rope("normal"), c(-0.10, 0.10))

    # Poisson / count / rate ratio
    expect_equal(get_family_default_rope("poisson"), c(-0.05, 0.05))
    expect_equal(get_family_default_rope("negbinomial"), c(-0.05, 0.05))
    expect_equal(get_family_default_rope("zip"), c(-0.05, 0.05))
    expect_equal(get_family_default_rope("zinb"), c(-0.05, 0.05))

    # Binomial / categorical / latent logistic
    expect_equal(get_family_default_rope("binomial"), c(-0.18, 0.18))
    expect_equal(get_family_default_rope("bernoulli"), c(-0.18, 0.18))
    expect_equal(get_family_default_rope("ordinal"), c(-0.18, 0.18))
    expect_equal(get_family_default_rope("multinomial"), c(-0.18, 0.18))

    # Fallback for unknown / custom
    expect_equal(get_family_default_rope("custom_dist"), c(-0.10, 0.10))
    expect_equal(get_family_default_rope(NULL), c(-0.10, 0.10))
    expect_equal(get_family_default_rope(NA), c(-0.10, 0.10))
})

test_that("resolve_rope_for_param prioritizes correctly", {
    # 1. Global numeric vector overrides all
    expect_equal(
        resolve_rope_for_param("beta_Y_X", resp = "Y", fam = "poisson", rope_spec = c(-0.2, 0.2)),
        c(-0.2, 0.2)
    )

    # 2. Exact parameter match in list
    spec_param <- list(beta_Y_X = c(-0.01, 0.01), Y = c(-0.05, 0.05), poisson = c(-0.08, 0.08))
    expect_equal(
        resolve_rope_for_param("beta_Y_X", resp = "Y", fam = "poisson", rope_spec = spec_param),
        c(-0.01, 0.01)
    )

    # 3. Response match in list
    spec_resp <- list(Y = c(-0.04, 0.04), poisson = c(-0.08, 0.08))
    expect_equal(
        resolve_rope_for_param("beta_Y_Z", resp = "Y", fam = "poisson", rope_spec = spec_resp),
        c(-0.04, 0.04)
    )

    # 4. Family match in list
    spec_fam <- list(poisson = c(-0.07, 0.07), default = c(-0.15, 0.15))
    expect_equal(
        resolve_rope_for_param("beta_Count_X", resp = "Count", fam = "poisson", rope_spec = spec_fam),
        c(-0.07, 0.07)
    )

    # 5. Default key in list
    spec_default <- list(gaussian = c(-0.05, 0.05), default = c(-0.12, 0.12))
    expect_equal(
        resolve_rope_for_param("beta_Count_X", resp = "Count", fam = "poisson", rope_spec = spec_default),
        c(-0.12, 0.12)
    )

    # 6. Fallback to principled family default
    expect_equal(
        resolve_rope_for_param("beta_Prob_X", resp = "Prob", fam = "binomial", rope_spec = NULL),
        c(-0.18, 0.18)
    )
})

test_that("explain_rope_spec produces informative explanations and notices", {
    # Empty
    expect_equal(explain_rope_spec(NULL)$lines, character(0))
    expect_equal(explain_rope_spec(data.frame())$lines, character(0))

    # Resolved df
    res_df <- data.frame(
        Response = c("mass", "counts", "binary"),
        Family   = c("gaussian", "poisson", "binomial"),
        r_min    = c(-0.10, -0.05, -0.18),
        r_max    = c(0.10, 0.05, 0.18),
        stringsAsFactors = FALSE
    )

    exp_res <- explain_rope_spec(res_df)
    expect_length(exp_res$lines, 3)
    expect_match(exp_res$lines[1], "Gaussian")
    expect_match(exp_res$lines[2], "Poisson")
    expect_match(exp_res$lines[3], "Binomial")

    # Data warning for unstandardized Gaussian
    df_raw <- data.frame(
        mass = c(10, 20, 30, 40, 50),
        counts = c(1, 2, 3, 4, 5),
        binary = c(0, 1, 0, 1, 1)
    )
    exp_warn <- explain_rope_spec(res_df, data = df_raw)
    expect_false(is.null(exp_warn$note))
    expect_match(exp_warn$note, "mass")

    # Shared specification across multiple responses (e.g. d-sep tests)
    df_shared <- data.frame(
        Response = c("BM", "DD", "LS", "NL", "RS"),
        Family   = rep("gaussian", 5),
        r_min    = rep(-0.10, 5),
        r_max    = rep(0.10, 5),
        stringsAsFactors = FALSE
    )
    exp_shared <- explain_rope_spec(df_shared)
    expect_length(exp_shared$lines, 1)
    expect_match(exp_shared$lines[1], "^Gaussian: \\[-0.10, 0.10\\]")
})
