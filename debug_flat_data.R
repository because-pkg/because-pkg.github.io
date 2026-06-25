devtools::load_all("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because")

set.seed(123)
data_list <- list(
    species = data.frame(Species = 1:50, Body_Mass_s = rnorm(50), Metabolic_Rate = rnorm(50), Thermal_Tol = rnorm(50)),
    site = data.frame(Site = 1:30, Elevation_s = rnorm(30), NDVI = rnorm(30), Flower_Cover = rnorm(30)),
    survey = data.frame(Survey = 1:90, Site = rep(1:30, each=3), Temperature = rnorm(90), Wind_Speed = rnorm(90)),
    obs = data.frame(Abundance = rpois(4500, 5), Survey = rep(1:90, each=50), Species = rep(1:50, times=90))
)

structures <- list(
    species = diag(50),
    site = diag(30)
)
names(structures$species) <- 1:50
names(structures$site) <- 1:30
class(structures$species) <- "matrix"
class(structures$site) <- "matrix"

eqs <- list(
    Abundance ~ Body_Mass_s + Thermal_Tol + (1 | Site) + (1 | Survey) + (1 | Species)
)

# We will patch `flatten_for_python` temporarily to print names
trace(because:::flatten_for_python, exit = quote(print(names(returnValue()))))

tryCatch({
    because(
         equations = eqs,
         data      = data_list,
         levels    = list(
             species = c("Body_Mass_s", "Metabolic_Rate", "Thermal_Tol", "Species"),
             site    = c("Elevation_s", "NDVI", "Flower_Cover", "U_Resource", "Site"),
             survey  = c("Temperature", "Wind_Speed", "Survey"),
             obs     = c("Abundance")),
         hierarchy    = "site > survey > obs; species > obs",
         link_vars    = list(site = "Site", survey = "Survey", species = "Species"),
         latent       = "U_Resource",
         family       = c(Abundance = "poisson"),
         structure    = structures,
         id_col       = "Species",
         dsep         = TRUE,
         engine     = "numpyro"
    )
}, error=function(e) print(e))
