library(cowplot)
library(ggplot2)

setwd("path/to/base-dir")

data <- data.frame(counts = scan("cell_counts.txt"))

p1 <- ggplot(data, aes(counts)) +
  geom_histogram(color = "white", fill = "orange") +
  scale_x_log10() +
  labs(x = "Cells per Sample", y = "Samples") +
  theme_cowplot()

ggsave(p1, file = "histogram.png", width = 2715, height = 1306, units = "px")
