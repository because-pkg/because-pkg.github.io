# Test Hierarchical Polynomials Debug

devtools::load_all(".")
library(ape)

set.seed(123)

# Level 1: Species (N=10)
N_Species <- 10
Traits_df <- data.frame(
    Species = paste0("sp", 1:N_Species),
    Trait = rnorm(N_Species)
)

# Level 2: Obs (N=50)
N_Obs <- 50
Y_df <- data.frame(
    Species = sample(Traits_df$Species, N_Obs, replace = TRUE),
    Y = rnorm(N_Obs),
    Age = runif(N_Obs, 1, 10)
)

print("--- Testing Manual Data Prep ---")
# 1. Manually setup inputs to prepare_hierarchical_jags_data
h_info <- list(
    data = list(Traits = Traits_df, Obs = Y_df),
    levels = list(Traits = c("Trait"), Obs = c("Y", "Age")),
    hierarchy = "Traits > Obs",
    link_vars = "Species"
)

# NOTE: We simulate the case where Polynomials are Auto-Generated
# BUT for now let's assume we need 'Age' to be passed so JAGS can compute Age^2
eq_vars <- c("Y", "Trait", "Age")

prep_res <- prepare_hierarchical_jags_data(h_info, eq_vars)
print("Data List names from Prep:")
print(names(prep_res$data_list))

if ("Age" %in% names(prep_res$data_list)) {
    print("SUCCESS: 'Age' found in data list")
} else {
    print("FAILURE: 'Age' NOT found in data list")
}

print("--- Running because() ---")
tryCatch(
    {
        res <- because(
            equations = list(Y ~ Trait + I(Age^2)),
            data = list(Traits = Traits_df, Obs = Y_df),
            quiet = FALSE
        )
        print("SUCCESS: Model compiled!")
    },
    error = function(e) {
        print(paste("FAILURE:", e$message))
    }
)
