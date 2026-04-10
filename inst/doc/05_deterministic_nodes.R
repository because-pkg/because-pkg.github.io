## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(because)
set.seed(123)

## ----ex1----------------------------------------------------------------------
N <- 100

# Predictors
Rain <- rnorm(N)
Temp <- rnorm(N)

# Interaction: Growth depends on Rain * Temp
Growth <- 0.5 * (Rain * Temp) + rnorm(N, sd = 0.1)

# remove some random values in Rain and Temp
Rain[sample(1:N, 10)] <- NA
Temp[sample(1:N, 15)] <- NA


data <- data.frame(Growth = Growth, Rain = Rain, Temp = Temp)

# Note the formula: Growth ~ Rain * Temp
fit_int <- because(
    equations = list(
        Growth ~ Rain * Temp
    ),
    data = data,
    n.iter = 1000,
    quiet = TRUE
)
# The model creates a deterministic node for the interaction

# Notice the parameter name "beta_Growth_Rain_x_Temp"
summary(fit_int)

## ----ex2----------------------------------------------------------------------
# Predictors
Age <- runif(N, 0, 10)

# Deterministic threshold: 1 if Age > 2, else 0
IsMature <- as.numeric(Age > 2)

# Outcome: Mating success depends on Maturity
Mating <- 2 * IsMature + rnorm(N, sd = 0.5)

data_logic <- data.frame(Mating = Mating, Age = Age)

# Use I() to wrap the logic
fit_logic <- because(
    equations = list(
        Mating ~ I(Age > 3)
    ),
    data = data_logic,
    n.iter = 1000,
    quiet = TRUE
)

# The model recovers the effect of the threshold
# Parameter: "beta_Mating_Age_gt_2" where "gt" means "greater than"
summary(fit_logic)

## ----ex3----------------------------------------------------------------------
# Define 4-level class: 1=Newborn, 2=Juv, 3=Sub, 4=Adult
# IMPORTANT: The class must be derived from Age so the model understands the structure.
Age <- round(runif(N, 0, 10), 0)

AgeClass <- 1 * (Age == 0) +
    2 * (Age > 0 & Age < 2) +
    3 * (Age >= 2 & Age < 5) +
    4 * (Age >= 5)

# Outcome: Social Rank increases with Life Stage
Rank <- 1.5 * AgeClass + rnorm(N, sd = 0.5)

data_multi <- data.frame(Rank = Rank, Age = Age, AgeClass = AgeClass)

# because() recovers the slope (~1.5) and confirms the causal path Age -> AgeClass -> Rank
fit_multi <- because(
    equations = list(
        # 1. Link AgeClass to Rank (Causal: Class -> Rank)
        Rank ~ AgeClass
    ),
    data = data_multi,
    n.iter = 1000,
    quiet = TRUE
)

# Check effect of AgeClass on Rank
summary(fit_multi)

## ----ex4----------------------------------------------------------------------
Mass <- runif(N, 1, 100) # Ensure positive for log
Metabolism <- 0.75 * log(Mass) + rnorm(N, sd = 0.1)

data_math <- data.frame(Metabolism = Metabolism, Mass = Mass)

fit_math <- because(
    equations = list(
        Metabolism ~ log(Mass)
    ),
    data = data_math,
    n.iter = 1000,
    quiet = TRUE
)

# Parameter: "beta_Metabolism_log_Mass"
summary(fit_math)

## -----------------------------------------------------------------------------
fit_multi_dsep <- because(
    equations = list(
        AgeClass ~ I(
            1 * (Age < 1) +
                2 * (Age >= 1 & Age < 2) +
                3 * (Age >= 2 & Age < 5) +
                4 * (Age >= 5)
        ),
        Rank ~ AgeClass
    ),
    data = data_multi,
    n.iter = 1000,
    dsep = TRUE,
    quiet = TRUE
)

summary(fit_multi_dsep)

