#!/usr/bin/env Rscript

library(argparser)
library(cowplot)
library(ggplot2)


p <- arg_parser("Plot Fig_03")

p <- add_argument(p, "--input", help = "path to input cell counts")
p <- add_argument(p, "--output", help = "output path for the plot")

argv <- parse_args(p)

input <- argv$input
output <- argv$output


data <- data.frame(counts = scan(input))

p1 <- ggplot(data, aes(counts)) +
  geom_histogram(color = "white", fill = "orange") +
  scale_x_log10() +
  labs(x = "Cells per Sample", y = "Samples") +
  theme_cowplot()

p1
ggsave(p1, file = output, width = 2715, height = 1306, units = "px")