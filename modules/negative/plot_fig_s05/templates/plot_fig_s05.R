#!/usr/bin/env Rscript
        
library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)

data.sim <- fread("${plotting_data_atlas_negative}")
data.real <- fread("${plotting_data_kang}")

data.sim <- melt(data.sim, id.vars = "pvals", variable.name = "method", value.name = "fps")
data.real <- melt(data.real, id.vars = "pvals", variable.name = "method", value.name = "fps")

data.sim <- na.omit(data.sim)
data.sim <- data.sim[method != "distinct"]

data.real <- na.omit(data.real)
data.real <- data.real[method != "distinct"]


color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation-test", "scvi"), 
                        color = c(1:6), 
                        hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
                        method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI"))

data.sim\$method <- factor(data.sim\$method, levels = color.code\$method)
data.real\$method <- factor(data.real\$method, levels = color.code\$method)


p1.sim <- ggplot(data.sim, aes(x = pvals, y = fps, color = method)) +
geom_line(linewidth= 0.8) +
coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 150)) +
labs(x = "p-value (except scVI)", y = "# False Positives", color = 'Method') +
scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
theme_cowplot(14) +
theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
background_grid(major = "xy",
                minor = "none",
                size.major = 0.2,
                size.minor = 0.2)


p2.sim <- ggplot(data.sim, aes(x = pvals, y = fps, color = method)) +
geom_line(linewidth= 0.8) +
labs(x = "p-value (except scVI)", y = "# False Positives", color = 'Method') +
scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
theme_cowplot(14) +
theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
background_grid(major = "xy",
                minor = "none",
                size.major = 0.2,
                size.minor = 0.2)


p1.real <- ggplot(data.real, aes(x = pvals, y = fps, color = method)) +
geom_line(linewidth= 0.8) +
coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 150)) +
labs(x = "p-value (except scVI)", y = "# False Positives", color = 'Method') +
scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
theme_cowplot(14) +
theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
background_grid(major = "xy",
                minor = "none",
                size.major = 0.2,
                size.minor = 0.2)


p2.real <- ggplot(data.real, aes(x = pvals, y = fps, color = method)) +
geom_line(linewidth= 0.8) +
labs(x = "p-value (except scVI)", y = "# False Positives", color = 'Method') +
scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
theme_cowplot(14) +
theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
background_grid(major = "xy",
                minor = "none",
                size.major = 0.2,
                size.minor = 0.2)

# Combine plots into one 2x2
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

ggsave(plot, filename = "Fig_S05.png", width = 3606, height = 4096, units = "px", dpi=400)