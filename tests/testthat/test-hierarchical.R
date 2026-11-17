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

test_that("because() preserves distinct multiscale structures rather than flattening to an ensemble", {
    testthat::skip_if_not_installed("rjags")
    testthat::skip_if_not_installed("ape")
    testthat::skip_if_not_installed("because.phybase")

    library(because.phybase)
    set.seed(123)
    n_sp <- 8
    n_site <- 8

    tree <- ape::rtree(n_sp)
    tree$tip.label <- paste0("Sp_", seq_len(n_sp))
    tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))

    coords <- cbind(runif(n_site), runif(n_site))
    rownames(coords) <- paste0("Site_", seq_len(n_site))
    v_spatial <- exp(-as.matrix(dist(coords)))

    d_sp <- data.frame(
        Species = tree$tip.label,
        Trait = rnorm(n_sp)
    )
    d_site <- data.frame(
        Site = rownames(coords),
        Env = rnorm(n_site)
    )
    d_obs <- expand.grid(Species = d_sp$Species, Site = d_site$Site, stringsAsFactors = FALSE)
    d_obs$obs_id <- seq_len(nrow(d_obs))
    d_obs$Y <- rnorm(nrow(d_obs))

    data_list <- list(species = d_sp, site = d_site, obs = d_obs)
    fit <- because(
        equations = list(Y ~ Trait + Env),
        data = data_list,
        levels = list(
            species = c("Trait", "Species"),
            site = c("Env", "Site"),
            obs = c("Y")
        ),
        hierarchy = "site > obs; species > obs",
        link_vars = list(site = "Site", species = "Species"),
        id_col = "Species",
        structure = list(phylo = tree, spatial = v_spatial),
        engine = "jags",
        n.chains = 2,
        n.iter = 100,
        n.burnin = 20,
        quiet = TRUE
    )

    expect_s3_class(fit, "because")
    expect_equal(names(fit$structures), c("phylo", "spatial"))
    expect_s3_class(fit$structures$phylo, "phylo")
    expect_true(is.matrix(fit$structures$spatial))
})
