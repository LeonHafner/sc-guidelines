#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)


meta_string <- "${meta_string}"
path_string <- "${path_string}"

parse_meta <- function(meta) {
  parts <- strsplit(meta, ", ")[[1]]
  scenario <- sub("scenario:", "", parts[1])
  run <- sub("run:", "", parts[2])
  return(list("scenario" = scenario, "run" = run))
}

# Preprocess meta_string and path_string
meta_components <- strsplit(gsub("\\\\[|\\\\]", "", meta_string), ";")[[1]]
path_components <- strsplit(path_string, ";")[[1]]

data <- lapply(1:length(meta_components), function(i) {
  meta_parsed <- parse_meta(meta_components[i])
  list("scenario" = meta_parsed\$scenario,
       "run" = meta_parsed\$run,
       "path" = path_components[i])
})
file.info <- do.call(rbind, lapply(data, as.data.frame))


# Load data
df <- data.table(method = character(), auc = numeric(), scenario = character())

for (scenario in c("atlas", "atlas-ub-conditions", "dataset", "dataset-ub-cells")) {
  
  files.auc <- file.info[file.info\$scenario == scenario,]\$path

  dfs <- list()
  for (file in files.auc) {
    dfs <- append(dfs, list(fread(file)))
  }
  
  df.per.scenario <- rbindlist(dfs)
  df.per.scenario\$scenario <- scenario
  
  df <- rbind(df, df.per.scenario)
}

df <- df[method != "distinct"]

color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation-test", "scvi"), 
                         color = c(1:6), 
                         hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"),
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI"))


# Atlas
df.atlas <- df[scenario == "atlas"]
df.atlas.ordered <- df.atlas[,.(median_auc = median(auc)), by = method][order(median_auc, decreasing = T)]

color.code.atlas <- merge(df.atlas.ordered, color.code, by = "method", sort = F)

df.atlas\$method <- factor(df.atlas\$method, levels = color.code.atlas\$method)

p.atlas <- ggplot(df.atlas, aes(x = method, y = auc, color = method)) +
  geom_boxplot() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(x = "Method", y = "AUPRC", color = "Method") +
  scale_color_okabe_ito(order = color.code.atlas\$color, labels = c(color.code.atlas\$method_legend)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "black", size = 10, hjust = 0.5)) +
  ggtitle("Atlas")


# Atlas-ub-conditions
df.atlas_ub_conditions <- df[scenario == "atlas-ub-conditions"]
df.atlas_ub_conditions.ordered <- df.atlas_ub_conditions[,.(median_auc = median(auc)), by = method][order(median_auc, decreasing = T)]

color.code.atlas_ub_conditions <- merge(df.atlas_ub_conditions.ordered, color.code, by = "method", sort = F)

df.atlas_ub_conditions\$method <- factor(df.atlas_ub_conditions\$method, levels = color.code.atlas_ub_conditions\$method)

p.atlas_ub_conditions <- ggplot(df.atlas_ub_conditions, aes(x = method, y = auc, color = method)) +
  geom_boxplot() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(x = "Method", y = "AUPRC", color = "Method") +
  scale_color_okabe_ito(order = color.code.atlas_ub_conditions\$color, labels = color.code.atlas_ub_conditions\$method_legend) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "black", size = 10, hjust = 0.5)) +
  ggtitle("Atlas with unbalanced conditions")

# Dataset
df.dataset <- df[scenario == "dataset"]
df.dataset.ordered <- df.dataset[,.(median_auc = median(auc)), by = method][order(median_auc, decreasing = T)]

color.code.dataset <- merge(df.dataset.ordered, color.code, by = "method", sort = F)

df.dataset\$method <- factor(df.dataset\$method, levels = color.code.dataset\$method)

p.dataset <- ggplot(df.dataset, aes(x = method, y = auc, color = method)) +
  geom_boxplot() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(x = "Method", y = "AUPRC", color = "Method") +
  scale_color_okabe_ito(order = color.code.dataset\$color, labels = color.code.dataset\$method_legend) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "black", size = 10, hjust = 0.5)) +
  ggtitle("Dataset")

# Dataset-ub-cells
df.dataset_ub_cells <- df[scenario == "dataset-ub-cells"]
df.dataset_ub_cells.ordered <- df.dataset_ub_cells[,.(median_auc = median(auc)), by = method][order(median_auc, decreasing = T)]

color.code.dataset_ub_cells <- merge(df.dataset_ub_cells.ordered, color.code, by = "method", sort = F)

df.dataset_ub_cells\$method <- factor(df.dataset_ub_cells\$method, levels = color.code.dataset_ub_cells\$method)

p.dataset_ub_cells <- ggplot(df.dataset_ub_cells, aes(x = method, y = auc, color = method)) +
  geom_boxplot() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(x = "Method", y = "AUPRC", color = "Method") +
  scale_color_okabe_ito(order = color.code.dataset_ub_cells\$color, labels = color.code.dataset_ub_cells\$method_legend) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(color = "black", size = 10, hjust = 0.5)) +
  ggtitle("Dataset with varying cell counts")

grid <- plot_grid(
  plot_grid(p.dataset + theme(legend.position = "none", plot.margin = unit(c(0.5, 0.1, 0, 0.1), "cm")),
            p.atlas + theme(legend.position = "none", axis.title.y = element_blank(), plot.margin = unit(c(0.5, 0.1, 0, 0.1), "cm")),
            p.dataset_ub_cells + theme(legend.position = "none", plot.margin = unit(c(0.5, 0.1, 0, 0.1), "cm")),
            p.atlas_ub_conditions + theme(legend.position = "none", axis.title.y = element_blank(), plot.margin = unit(c(0.5, 0.1, 0, 0.1), "cm")),
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

ggsave(plot = grid, filename = "Fig_05.png", width = 1803, height = 2048, units = "px")

df.atlas.ordered[, scenario := "atlas"]
df.atlas_ub_conditions.ordered[, scenario := "atlas_ub_conditions"]
df.dataset.ordered[, scenario := "dataset"]
df.dataset_ub_cells.ordered[, scenario := "dataset_ub_cells"]

fwrite(rbindlist(list(df.atlas.ordered, df.atlas_ub_conditions.ordered, df.dataset.ordered, df.dataset_ub_cells.ordered)), file = "median_auc.tsv", sep = "\\t")
