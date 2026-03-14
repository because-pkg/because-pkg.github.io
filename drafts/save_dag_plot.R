options(digits = 3)
devtools::load_all("../because")
devtools::load_all("../because.phybase")
library(ggplot2)

# Load the fitted simulation
fit <- readRDS("simulated_dual_scale.rds")

# Generate the DAG plot (Returns a ggplot object)
my_dag <- plot_dag(fit)

# Export the network to PNG
ggsave(
    "drafts/simulated_dag.png",
    plot = my_dag,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
)
cat("DAG plot saved to drafts/simulated_dag.png\n")
