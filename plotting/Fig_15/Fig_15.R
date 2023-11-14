library(data.table)
library(cowplot)
library(ggplot2)
library(ggokabeito)


setwd('path/to/base-dir')

path.prc <- file.path(c('atlas', 'atlas-ub-conditions', 'dataset', 'dataset-ub-cells'),
                      'data/13_precision-recall')
# Use arbitrary simulation run
files <- 'prc_distinct_5.tsv'
prc <- data.table(precision = numeric(), recall = numeric(), method = character(), scenario = factor())

# Read precision-recall data for each scenario and method
for (path in path.prc) {
  for (file in files) {
    data.method <- fread(file.path(path, file))
    data.method$method <- strsplit(strsplit(file, split = '.', fixed = T)[[1]][1], split = '_', fixed = T)[[1]][[2]]
    data.method$scenario <- strsplit(path, split = "/", fixed = T)[[1]][8]
    prc <- rbind(prc, data.method)
  }
}

# Read AUPRC data for each scenario
auc <- data.table(method = character(), auc = numeric(), scenario = character())
for (scenario in c("atlas", "atlas-ub-conditions", "dataset", "dataset-ub-cells")) {
  
  path <- file.path(base.path, scenario, 'data/13_precision-recall')
  files.auc <- list.files(path)
  files.auc <- files.auc[startsWith(files.auc, "auc")]
  
  aucs <- list()
  for (file in files.auc) {
    aucs <- append(aucs, list(fread(file.path(path, file))))
  }
  
  auc.per.scenario <- rbindlist(aucs)
  auc.per.scenario$scenario <- scenario
  
  auc <- rbind(auc, auc.per.scenario)
}

# Filter methods to distinct
auc <- auc[method == "distinct"]

# Set color coding for scenarios
color.code <- data.table(scenario = c("atlas", "atlas-ub-conditions", "dataset", "dataset-ub-cells"), 
                         color = c(1:4), 
                         scenario_legend = c("Atlas", "Atlas (unbalanced conditions)", "Dataset", "Dataset (unbalanced cells)"))

# Order by AUC (decreasing)
auc.ordered <- auc[,.(median_auc = median(auc)), by = "scenario"][order(median_auc, decreasing = T)]
color.code.auc <- merge(auc.ordered, color.code, by = "scenario", sort = F)

# Factorize to set legend order
prc$scenario <- factor(prc$scenario, levels = color.code.auc$scenario)
auc$scenario <- factor(auc$scenario, levels = color.code.auc$scenario)

# Plot PRC data
plot.prc <- ggplot(prc, aes(x = recall, y = precision, color = scenario)) +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
  geom_line() +
  ylim(0, 1) +
  scale_color_okabe_ito(order = color.code.auc$color, labels = color.code.auc$scenario_legend) +
  labs(x = "Recall", y = "Precision", color = "Scenario") +
  theme_cowplot()


# Plot AUPRC data
plot.auc <- ggplot(auc, aes(x = scenario, y = auc, color = scenario)) +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") + 
  geom_boxplot() +
  ylim(0, 1) +
  labs(color = "Scenario", y = "AUPRC") +
  scale_color_okabe_ito(order = color.code.auc$color, labels = color.code.auc$scenario_legend) +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Combine plots into 1x2 grid
p <- plot_grid(
              plot_grid(plot.prc + theme(legend.position = "none"),
                        plot.auc + theme(legend.position = "none"),
                        labels = "AUTO",
                        label_size = 20,
                        ncol = 2,
                        label_x = 0,
                        label_y = 1,
                        hjust = -0.1,
                        vjust = 1.1,
                        align = "hv",
                        rel_widths = c(1, 1.1)),
              get_legend(plot.prc +
                           theme(legend.text = element_text(size = 10),
                                 legend.title = element_text(size = 14),
                                 legend.spacing.x = unit(1.0, "cm"),
                                 legend.title.align = 0.5,
                                 legend.position=c(0.2, 0.6)) +
                           guides(color = guide_legend(ncol = 2))),
              ncol = 1,
              rel_heights = c(1, 0.2)
              )

ggsave(p, width = 3606, height = 1950, units = "px", dpi = 400, filename = "Fig_15.png")




