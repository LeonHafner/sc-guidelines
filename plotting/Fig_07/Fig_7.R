library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)

setwd('path/to/base-dir')

df <- data.table(method = character(), auc = numeric(), scenario = character())

for (scenario in c("atlas", "atlas-ub-conditions", "dataset", "dataset-ub-cells")) {
  
  path <- file.path(scenario, 'data/13_precision-recall')
  files.auc <- list.files(path)
  files.auc <- files.auc[startsWith(files.auc, "auc")]
  
  dfs <- list()
  for (file in files.auc) {
    dfs <- append(dfs, list(fread(file.path(path, file))))
  }
  
  df.per.scenario <- rbindlist(dfs)
  df.per.scenario$scenario <- scenario
  
  df <- rbind(df, df.per.scenario)
}

df <- df[method != "distinct"]


color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation", "scvi"), 
                         color = c(1:6), 
                         hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI"))


# Plotting the atlas scenario
df.atlas <- df[scenario == "atlas"]
df.atlas.ordered <- df.atlas[,.(median_auc = median(auc)), by = method][order(median_auc, decreasing = T)]

color.code.atlas <- merge(df.atlas.ordered, color.code, by = "method", sort = F)
color.code.atlas

df.atlas$method <- factor(df.atlas$method, levels = color.code.atlas$method)

p.atlas <- ggplot(df.atlas, aes(x = method, y = auc, color = method)) +
  geom_boxplot() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(x = "Method", y = "AUPRC", color = "Method") +
  scale_color_okabe_ito(order = color.code.atlas$color, labels = c(color.code.atlas$method_legend)) +
  #scale_x_discrete(labels=color.code.atlas$method_legend) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())


# Plotting the atlas scenario with unbalanced conditions
df.atlas_ub_conditions <- df[scenario == "atlas-ub-conditions"]
df.atlas_ub_conditions.ordered <- df.atlas_ub_conditions[,.(median_auc = median(auc)), by = method][order(median_auc, decreasing = T)]

color.code.atlas_ub_conditions <- merge(df.atlas_ub_conditions.ordered, color.code, by = "method", sort = F)
color.code.atlas_ub_conditions

df.atlas_ub_conditions$method <- factor(df.atlas_ub_conditions$method, levels = color.code.atlas_ub_conditions$method)

p.atlas_ub_conditions <- ggplot(df.atlas_ub_conditions, aes(x = method, y = auc, color = method)) +
  geom_boxplot() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(x = "Method", y = "AUPRC", color = "Method") +
  scale_color_okabe_ito(order = color.code.atlas_ub_conditions$color, labels = color.code.atlas_ub_conditions$method_legend) +
  #scale_x_discrete(labels=color.code.atlas_ub_conditions$method_legend) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p.atlas_ub_conditions


# Plotting the dataset scenario
df.dataset <- df[scenario == "dataset"]
df.dataset.ordered <- df.dataset[,.(median_auc = median(auc)), by = method][order(median_auc, decreasing = T)]

color.code.dataset <- merge(df.dataset.ordered, color.code, by = "method", sort = F)
color.code.dataset

df.dataset$method <- factor(df.dataset$method, levels = color.code.dataset$method)

p.dataset <- ggplot(df.dataset, aes(x = method, y = auc, color = method)) +
  geom_boxplot() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(x = "Method", y = "AUPRC", color = "Method") +
  scale_color_okabe_ito(order = color.code.dataset$color, labels = color.code.dataset$method_legend) +
  #scale_x_discrete(labels=color.code.dataset$method_legend) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())


# Plotting the dataset scenario with unbalanced cell numbers
df.dataset_ub_cells <- df[scenario == "dataset-ub-cells"]
df.dataset_ub_cells.ordered <- df.dataset_ub_cells[,.(median_auc = median(auc)), by = method][order(median_auc, decreasing = T)]

color.code.dataset_ub_cells <- merge(df.dataset_ub_cells.ordered, color.code, by = "method", sort = F)
color.code.dataset_ub_cells

df.dataset_ub_cells$method <- factor(df.dataset_ub_cells$method, levels = color.code.dataset_ub_cells$method)

p.dataset_ub_cells <- ggplot(df.dataset_ub_cells, aes(x = method, y = auc, color = method)) +
  geom_boxplot() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(x = "Method", y = "AUPRC", color = "Method") +
  scale_color_okabe_ito(order = color.code.dataset_ub_cells$color, labels = color.code.dataset_ub_cells$method_legend) +
  #scale_x_discrete(labels=color.code.dataset_ub_cells$method_legend) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())
p.dataset_ub_cells


# Combining the plots into one 2x2 grid
grid <- plot_grid(
          plot_grid(p.atlas + theme(legend.position = "none"),
                    p.atlas_ub_conditions + theme(axis.title.y = element_blank(), legend.position = "none"),
                    p.dataset + theme(legend.position = "none"),
                    p.dataset_ub_cells + theme(axis.title.y = element_blank(), legend.position = "none"),
                    labels = "AUTO",
                    ncol = 2,
                    align = "hv",
                    axis = "tb",
                    rel_widths = c(1, 1.07, 1, 1.07)),
          get_legend(p.atlas +
                       theme(legend.text = element_text(size = 11),
                             legend.title = element_text(size = 15),
                             legend.spacing.x = unit(0.5, "cm"),
                             legend.title.align = 0.5,
                             legend.position=c(0.2, 0.6)) +
                       guides(color = guide_legend(ncol = 3))),
        nrow = 2,
        rel_heights = c(1, 0.3),
        vjust = -10)

ggsave(grid, "Fig_7.png", width = 5409, height = 6144, units = "px")
