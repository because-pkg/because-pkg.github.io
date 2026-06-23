library(ggplot2)
library(ggridges)

set.seed(1)
plot_data <- data.frame(
  Label = rep(c("A", "B", "C"), each = 1000),
  Value = c(rnorm(1000, 0, 0.05), rnorm(1000, 0.1, 0.1), rnorm(1000, -0.2, 0.05))
)

rope <- c(-0.1, 0.1)

p <- ggplot(plot_data, aes(x = Value, y = Label)) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    aes(fill = after_stat(x >= rope[1] & x <= rope[2])),
    color = "black",
    scale = 0.95,
    rel_min_height = 0.005
  ) +
  scale_fill_manual(
    name = "Region",
    values = c("TRUE" = "#89C5C5", "FALSE" = "#F8766D")
  ) +
  geom_vline(xintercept = rope, linetype = "dashed", color = "darkgray")

ggsave("test.png", p, width=6, height=4)
