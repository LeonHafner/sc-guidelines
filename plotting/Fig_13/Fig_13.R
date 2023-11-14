library(ggplot2)
library(data.table)
library(pracma)
library(cowplot)
library(tidyr)
library(ggokabeito)

setwd("path/to/base-dir")

path.prc.sim <- 'sim-data/data/13_negative-control'
path.prc.real <- 'real-data/data/13_negative-control'

output.sim <- 'sim-data/data/14_plots'
output.real <- 'real-data/data/14_plots'

methods <- c('mast', 'distinct', 'deseq2', 'permutation', 'hierarchical-bootstrapping', 'scvi', 'dream')

files <- paste0('fps_', methods, '.tsv')

# Processing for simulated data
cutoffs.sim <- data.table(cutoff = numeric(), false_positives = numeric(), method = character())
for (file in files) {
  method <- fread(file.path(path.prc.sim, file))
  method$method <- strsplit(strsplit(file, split = '.', fixed = T)[[1]][1], split = '_', fixed = T)[[1]][[2]]
  cutoffs.sim <- rbind(cutoffs.sim, method)
}
cutoffs.sim <- na.omit(cutoffs.sim)
cutoffs.sim <- cutoffs.sim[method != "distinct"]


# Processing for real data
cutoffs.real <- data.table(cutoff = numeric(), false_positives = numeric(), method = character())
for (file in files) {
  method <- fread(file.path(path.prc.real, file))
  method$method <- strsplit(strsplit(file, split = '.', fixed = T)[[1]][1], split = '_', fixed = T)[[1]][[2]]
  cutoffs.real <- rbind(cutoffs.real, method)
}
cutoffs.real <- na.omit(cutoffs.real)
cutoffs.real <- cutoffs.real[method != "distinct"]


# Color coding of methods
color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation", "scvi"), 
                         color = c(1:6), 
                         hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI"))
 

# Factorizing the method column for legend order
cutoffs.sim$method <- factor(cutoffs.sim$method, levels = color.code$method)
cutoffs.real$method <- factor(cutoffs.real$method, levels = color.code$method)


p1.sim <- ggplot(cutoffs.sim, aes(x = cutoff, y = false_positives, color = method)) +
  geom_line(linewidth= 0.8) +
  xlim(0, 0.05) +
  ylim(0, 150) +
  labs(x = "p-value cutoff", y = "# False Positives", color = 'Method') +
  scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
  theme_cowplot(14) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  background_grid(major = "xy",
                minor = "none",
                size.major = 0.2,
                size.minor = 0.2)


p2.sim <- ggplot(cutoffs.sim, aes(x = cutoff, y = false_positives, color = method)) +
  geom_line(linewidth= 0.8) +
  labs(x = "p-value cutoff", y = "# False Positives", color = 'Method') +
  scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
  theme_cowplot(14) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  background_grid(major = "xy",
                  minor = "none",
                  size.major = 0.2,
                  size.minor = 0.2)


p1.real <- ggplot(cutoffs.real, aes(x = cutoff, y = false_positives, color = method)) +
  geom_line(linewidth= 0.8) +
  xlim(0, 0.1) +
  ylim(0, 150) +
  labs(x = "p-value cutoff", y = "# False Positives", color = 'Method') +
  scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
  theme_cowplot(14) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  background_grid(major = "xy",
                  minor = "none",
                  size.major = 0.2,
                  size.minor = 0.2)


p2.real <- ggplot(cutoffs.real, aes(x = cutoff, y = false_positives, color = method)) +
  geom_line(linewidth= 0.8) +
  labs(x = "p-value cutoff", y = "# False Positives", color = 'Method') +
  scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
  theme_cowplot(14) +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  background_grid(major = "xy",
                  minor = "none",
                  size.major = 0.2,
                  size.minor = 0.2)


# Combine plots into 2x2 grid
combined.plot <- plot_grid(p1.sim + theme(legend.position = "none",
                                          axis.title.x = element_blank()),
                           p2.sim + theme(legend.position = "none",
                                          axis.title.x = element_blank(),
                                          axis.title.y = element_blank()),
                           p1.real + theme(legend.position = "none"),
                           p2.real + theme(legend.position = "none",
                                           axis.title.y = element_blank()),
                           labels = "AUTO",
                           label_size = 20,
                           ncol = 2,
                           label_x = 0,
                           label_y = 1,
                           hjust = -0.1,
                           vjust = 1.1,
                           align = "hv")

# Add shared legend
plot <- plot_grid(combined.plot,
                  get_legend(p1.sim +
                               theme(legend.title = element_text(size = 15),
                                     legend.text = element_text(size = 11),
                                     legend.spacing.x = unit(0.5, "cm"),
                                     legend.title.align = 0.5,
                                     legend.position=c(0.3, 0.6)) +
                               guides(color = guide_legend(ncol = 3))),
                  ncol = 1,
                  rel_heights = c(4, 1))

ggsave(plot, width = 3606, height = 4096, units = "px", dpi = 400, filename = "Fig_13.png")
