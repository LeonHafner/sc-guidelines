library(data.table)
library(ggplot2)
library(ggokabeito)
library(cowplot)

setwd("path/to/base-dir")

# Load data
data.sim <- fread("sim-data/data/13_negative-control/plotting_data.tsv")
data.sim <- melt(data.sim, id.vars = "pvals", variable.name = "method", value.name = "fpr")

data.real <- fread("real-data/data/13_negative-control/plotting_data.tsv")
data.real <- melt(data.real, id.vars = "pvals", variable.name = "method", value.name = "fpr")


# Color coding of the methods
color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation", "scvi", "distinct"), 
                         color = c(1:7), 
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "distinct"))


data.sim$method <- factor(data.sim$method, levels = color.code$method)
data.sim <- na.omit(data.sim)
data.sim <- data.sim[method != "distinct"]

data.real$method <- factor(data.real$method, levels = color.code$method)
data.real <- na.omit(data.real)
data.real <- data.real[method != "distinct"]


# Plotting simulated data
p.sim <- ggplot(data.sim, aes(x = pvals, y = fpr, color = method)) +
  geom_line(linewidth = 1) +
  scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
  labs(x = "p-value (adjusted with BH)", y = "FPR", color = "Method") +
  geom_abline(linetype = "dashed", color = "gray") +
  theme_cowplot()


# Plotting real data
p.real <- ggplot(data.real, aes(x = pvals, y = fpr, color = method)) +
  geom_line(linewidth = 1) +
  scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
  labs(x = "p-value (adjusted with BH)", y = "FPR", color = "Method") +
  geom_abline(linetype = "dashed", color = "gray") +
  theme_cowplot()

# Combine plots into one panel
double.plot <- plot_grid(p.sim + theme(legend.position = "none",
                                       axis.title.x = element_blank()),
                         p.real + theme(legend.position = "none"),
                         labels = "AUTO",
                         label_size = 20,
                         ncol = 1,
                         label_x = 0, label_y = 1,
                         hjust = -0.1, vjust = 1.1,
                         align = "v")

legend.double.plot <- get_legend(p.sim + theme(legend.margin = margin(0, 0, 20, 20)))
p.combined <- plot_grid(double.plot, legend.double.plot, rel_widths = c(2.5, 1))

ggsave(p.combined, width = 2370, height = 3194, units = "px", dpi = 300, filename = "Fig_12.png")
