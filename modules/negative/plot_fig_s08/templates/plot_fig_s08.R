#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggokabeito)
library(cowplot)

data.sim <- fread("${plotting_data_atlas_negative}")
data.real <- fread("${plotting_data_kang}")

data.sim <- melt(data.sim, id.vars = "pvals", variable.name = "method", value.name = "fpr")
data.real <- melt(data.real, id.vars = "pvals", variable.name = "method", value.name = "fpr")

color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation-test", "scvi", "distinct", "ttest"), 
                        color = c(1:8), 
                        method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "distinct", "t-test"))


# Factorize methods to get the same color across all plots
data.sim\$method <- factor(data.sim\$method, levels = color.code\$method)
data.sim <- na.omit(data.sim)
data.sim <- data.sim[method != "distinct"]

data.real\$method <- factor(data.real\$method, levels = color.code\$method)
data.real <- na.omit(data.real)
data.real <- data.real[method != "distinct"]


p.sim <- ggplot(data.sim, aes(x = pvals, y = fpr, color = method)) +
geom_line(linewidth = 1) +
scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
labs(x = "p-value (except scVI)", y = "FPR", color = "Method") +
geom_abline(linetype = "dashed", color = "gray") +
theme_cowplot()

p.real <- ggplot(data.real, aes(x = pvals, y = fpr, color = method)) +
geom_line(linewidth = 1) +
scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
labs(x = "p-value (except scVI)", y = "FPR", color = "Method") +
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

ggsave(p.combined, filename = "Fig_S08.png", width = 2370, height = 3194, units = "px")