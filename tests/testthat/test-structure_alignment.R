test_that("collect_and_filter_model_data preserves structure alignment columns", {
    # Mock data with 6 columns, where 'Species' matches tree tips
    df <- data.frame(
        Species = c("Sp1", "Sp2", "Sp3", "Sp4"),
        LS = rnorm(4),
        DD = rnorm(4),
        BM = rnorm(4),
        NL = rnorm(4),
        Unrelated = rnorm(4),
        stringsAsFactors = FALSE
    )

    # Mock structure (phylo-like or matrix)
    mock_tree <- structure(
        list(tip.label = c("Sp1", "Sp2", "Sp3", "Sp4")),
        class = "phylo"
    )

    # Test equation with only 4 columns (LS, DD, BM, NL)
    eqs <- list(LS ~ DD + BM + NL)

    filtered_res <- collect_and_filter_model_data(
        data = df,
        equations = eqs,
        random = NULL,
        id_col = NULL,
        link_vars = NULL,
        family = NULL,
        variability = NULL,
        structure = mock_tree,
        quiet = TRUE
    )

    # 'Species' must be retained for alignment even though it is not in eqs
    expect_true("Species" %in% names(filtered_res$data))
    # Unrelated column must be dropped
    expect_false("Unrelated" %in% names(filtered_res$data))
    expect_equal(sort(names(filtered_res$data)), sort(c("Species", "LS", "DD", "BM", "NL")))
})

test_that("get_order_labels_hook extracts tip labels from phylo objects", {
    mock_tree <- structure(
        list(tip.label = c("A", "B", "C")),
        class = "phylo"
    )
    expect_equal(get_order_labels_hook(mock_tree), c("A", "B", "C"))
})

test_that("preprocess_categorical_vars sanitizes column names with dots for dummy variables", {
    df <- data.frame(
        Y = rnorm(5),
        Sample.size = c("large", "medium", "small", "tiny", "large"),
        stringsAsFactors = FALSE
    )
    res <- preprocess_categorical_vars(df, target_vars = "Sample.size", dummy_vars = "Sample.size", quiet = TRUE)
    cat_vars <- attr(res, "categorical_vars")
    expect_true("Sample.size" %in% names(cat_vars))
    # Dummies must use clean name with underscore
    expect_equal(cat_vars[["Sample.size"]]$dummies, c("Sample_size_medium", "Sample_size_small", "Sample_size_tiny"))
    expect_true(all(c("Sample_size_medium", "Sample_size_small", "Sample_size_tiny") %in% names(res)))

    # Expansion in equations must use sanitized dummy names
    eqs <- list(Y ~ Sample.size)
    expanded <- expand_categorical_equations(eqs, res, quiet = TRUE)
    exp_terms <- attr(terms(expanded$equations[[1]]), "term.labels")
    expect_equal(exp_terms, c("Sample_size_medium", "Sample_size_small", "Sample_size_tiny"))
})
