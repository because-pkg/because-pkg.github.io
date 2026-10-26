test_that("because() handles complex hierarchical/multi-level data across scales", {
    set.seed(42)
    
    # 1. Macro-level (Year/Ecological level)
    n_years <- 10
    df_year <- data.frame(
        YearID = 1:n_years,
        Climate = rnorm(n_years)
    )
    
    # Yearly effect: 2.0 * Climate + noise
    YearEffect <- 2.0 * df_year$Climate + rnorm(n_years, 0, 0.5)
    
    # 2. Micro-level (Individual level with repeated measures)
    # 5 individuals per year
    n_inds_per_year <- 5
    n_total <- n_years * n_inds_per_year
    
    df_ind <- data.frame(
        IndID = 1:n_total,
        YearID = rep(1:n_years, each = n_inds_per_year),
        Trait = rnorm(n_total)
    )
    
    # Trait affects Response, but Response also affected by YearEffect
    # Response = 1.5 * Trait + YearEffect + noise
    df_ind$Response <- 1.5 * df_ind$Trait + YearEffect[df_ind$YearID] + rnorm(n_total, 0, 0.5)
    
    # Pass as a list of dataframes
    data_list <- list(
        Year = df_year,
        Individual = df_ind
    )
    
    # Run the model
    # Climate -> Response (mediated implicitly via YearEffect, or directly mapped)
    # Trait -> Response
    fit <- because(
        list(Response ~ Trait + Climate),
        data = data_list,
        engine = "jags",
        n.iter = 100,
        n.burnin = 10,
        quiet = TRUE
    )
    
    expect_s3_class(fit, "because")
    
    # Verify the model compiled properly with the cross-level mapping
    expect_true("beta_Response_Climate" %in% fit$monitor)
    expect_true("beta_Response_Trait" %in% fit$monitor)
})
