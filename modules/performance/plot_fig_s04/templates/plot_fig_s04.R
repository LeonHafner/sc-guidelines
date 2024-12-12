#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(data.table)
})


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

dataframe <- data.table(method = character(), auc = numeric(), scenario = character())

for (scenario in unique(file.info\$scenario)) {
  
  files.auc <- file.info[file.info\$scenario == scenario,]\$path
  
  dfs <- list()
  for (file in files.auc) {
    dfs <- append(dfs, list(fread(file)))
  }
  
  df.per.scenario <- rbindlist(dfs)
  df.per.scenario\$scenario <- scenario
  
  dataframe <- rbind(dataframe, df.per.scenario)
}

# Introduce less_de column to data frame
dataframe\$less_de <- grepl("-less-de", dataframe\$scenario)

# Remove less-de information stored in the scenario column
dataframe\$scenario <- gsub("-less-de", "", dataframe\$scenario)

# Factorize scenario for right order of the plots
dataframe\$scenario <- factor(dataframe\$scenario, levels = c("dataset", "atlas", "dataset-ub-cells", "atlas-ub-conditions"))

# Custom labeling function
custom_labels <- function(variable, value) {
  # Define custom labels based on 'category'
  labels <- c(
    dataset = "Dataset",
    atlas = "Atlas",
    dataset_ub_cells = "Dataset with varying cell numbers",
    atlas_ub_conditions = "Atlas with unbalanced conditions"
  )
  return(labels[value])
}


p <- ggplot(dataframe, aes(x = method, y = auc, fill = less_de)) +
  geom_boxplot() +
  ylim(0, 1) +
  scale_fill_manual(values = c("#FF7F00", "#1F78B4"), labels = c("5%", "0.5%")) +
  scale_x_discrete(labels = c("DESeq2", "distinct", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "t-test")) +
  labs(x = "Method", y = "AUPRC", fill = "Percentage\nof DE genes") +
  facet_wrap(~ scenario, labeller = labeller(scenario = custom_labels)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = p, width = 3606, height = 3900, units = "px", dpi = 400, filename = "Fig_S04.png")
