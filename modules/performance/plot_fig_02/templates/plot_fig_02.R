#!/usr/bin/env Rscript

library(cowplot)
library(ggplot2)
library(png)
library(grid)

# Download Fig02A
png_url = "https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_02/Fig_2.png"
drawio_url = "https://raw.githubusercontent.com/LeonHafner/sc-guidelines/2ceb8da423dba514f503959fd1a56643e389eaf7/plotting/Fig_02/Fig_2.drawio"

download.file(png_url, destfile = "Fig_02A.png")
download.file(drawio_url, destfile = "Fig_02A.drawio")

input <- "${cell_counts}"

img <- readPNG("Fig_02A.png")
p1 <- rasterGrob(img, interpolate = T)

data <- data.frame(counts = scan(input))

p2 <- ggplot(data, aes(counts)) +
  geom_histogram(color = "white", fill = "orange") +
  scale_x_log10() +
  labs(x = "Cells per Sample", y = "Samples") +
  theme_cowplot()

plot_combined <- plot_grid(p1, p2,
                           ncol=1,
                           labels = "AUTO",
                           label_size = 25,
                           label_x = 0,
                           label_y = 1,
                           hjust = -0.1,
                           vjust = 1.1,
                           rel_heights=c(1.5, 1)
                           )

ggsave(plot_combined, file = "Fig_02.png", width = 2715, height = 3700, units = "px")