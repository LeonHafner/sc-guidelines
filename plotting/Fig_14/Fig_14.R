library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)


# Load data for simulations with 0.5% DE
base.path <- 'path/to/base-dir/lessDE'
df <- data.table(method = character(), auc = numeric(), scenario = character(), less_de = logical())
for (scenario in c("atlas", "atlas-ub-conditions", "dataset", "dataset-ub-cells")) {
  
  path <- file.path(base.path, scenario, 'data/13_precision-recall')
  files.auc <- list.files(path)
  files.auc <- files.auc[startsWith(files.auc, "auc")]
  
  dfs <- list()
  for (file in files.auc) {
    dfs <- append(dfs, list(fread(file.path(path, file))))
  }
  
  df.per.scenario <- rbindlist(dfs)
  df.per.scenario$scenario <- scenario
  df.per.scenario$less_de <- TRUE
  
  df <- rbind(df, df.per.scenario)
}


# Load data for simulations with 5% DE
base.path <- 'path/to/base-dir/moreDE'
for (scenario in c("atlas", "atlas-ub-conditions", "dataset", "dataset-ub-cells")) {
  
  path <- file.path(base.path, scenario, 'data/13_precision-recall')
  files.auc <- list.files(path)
  files.auc <- files.auc[startsWith(files.auc, "auc")]
  
  dfs <- list()
  for (file in files.auc) {
    dfs <- append(dfs, list(fread(file.path(path, file))))
  }
  
  df.per.scenario <- rbindlist(dfs)
  df.per.scenario$scenario <- scenario
  df.per.scenario$less_de <- FALSE
  
  df <- rbind(df, df.per.scenario)
}


# Custom labeling function
custom_labels <- function(variable, value) {
  # Define custom labels based on 'category'
  labels <- c(
    atlas = "Atlas",
    atlas_ub_conditions = "Atlas (unbalanced conditions)",
    dataset = "Dataset",
    dataset_ub_cells = "Dataset (unbalanced cells)"
  )
  return(labels[value])
}

p <- ggplot(df, aes(x = method, y = auc, fill = less_de)) +
  geom_boxplot() +
  ylim(0, 1) +
  scale_fill_manual(values = c("#FF7F00", "#1F78B4"), labels = c("5%", "0.5%")) +
  scale_x_discrete(labels = c("DESeq2", "distinct", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI")) +
  labs(x = "Method", y = "AUPRC", fill = "Percentage\nof DE genes") +
  facet_wrap(~ scenario, labeller = labeller(scenario = custom_labels)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(p, width = 3606, height = 3900, units = "px", dpi = 400, filename = "Fig_14.png")








