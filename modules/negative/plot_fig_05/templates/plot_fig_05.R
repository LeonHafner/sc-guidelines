#!/usr/bin/env Rscript
        
library(data.table)
library(ggplot2)
library(ggokabeito)
library(cowplot)

data.sim <- fread("${plotting_data}")
data.sim[1, 2:ncol(data.sim) := 0]
data.sim <- melt(data.sim, id.vars = "pvals", variable.name = "method", value.name = "fpr")

color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation-test", "scvi", "distinct", "ttest"), 
                         color = c(1, 2, 3, 4, 5, 6, 7, 9),
                         hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"),
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "distinct", "t-test"))


data.sim\$method <- factor(data.sim\$method, levels = color.code\$method)
data.sim <- na.omit(data.sim)
data.sim <- data.sim[method != "distinct"]

color.code <- color.code[method != "distinct"]

p.sim <- ggplot(data.sim, aes(x = pvals, y = fpr, color = method)) +
geom_line(linewidth = 1) +
scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
labs(x = "p-value (except scVI)", y = "FPR", color = "Method") +
geom_abline(linetype = "dashed", color = "gray") +
theme_cowplot() +
coord_cartesian(xlim = c(0, 0.15), ylim = c(0, 0.15))

ggsave(p.sim, width = 2650, height = 2000, units = "px", dpi = 300, filename = "Fig_05.png")