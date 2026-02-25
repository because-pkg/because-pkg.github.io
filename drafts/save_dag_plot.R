options(digits = 3)
devtools::load_all("../because")
devtools::load_all("../because.phybase")
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Load the fitted simulation
fit <- readRDS("simulated_dual_scale.rds")

# Generate the DAG plot natively
my_dag <- plot_dag(fit)

# Export the network to an SVG and then to PNG
temp_html <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(my_dag, file = temp_html, selfcontained = TRUE)
webshot2::webshot(
    temp_html,
    file = "drafts/simulated_dag.png",
    vwidth = 800,
    hwidth = 600,
    zoom = 2
)
cat("DAG plot saved to drafts/simulated_dag.png\n")
