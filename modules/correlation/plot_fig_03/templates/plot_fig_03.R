#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)


correlation_string <- "${correlation_string}"

print(correlation_string)

paths <- strsplit(correlation_string, ";")[[1]]

data <- rbindlist(lapply(paths, fread))

data[, run := factor(run)]

sample_size <- ceiling(nrow(data) * 0.01)
sampled_data <- data[sample(nrow(data), sample_size), ]

color.intra <- '#ff7f00'
color.inter <- '#1f78b4'
    
p1 <- ggplot(sampled_data, aes(x = run, y = correlation, fill = type)) +
  geom_boxplot() +
  labs(x = 'Simulation', y = 'Correlation (spearman)', fill = "Type of\ncorrelation") +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12)) +
  scale_x_discrete(labels = paste('Sim', 1:10, sep = '')) +
  scale_fill_manual(values = c(color.inter, color.intra)) +
  theme_cowplot(14)

ggsave(p1, file = "Fig_03.png", width = 2480, height = 1351, units = 'px')