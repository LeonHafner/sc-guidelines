#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)


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

for (scenario in unique(file.info\$scenario)) {
  
  files.auc <- file.info[file.info\$scenario == scenario,]\$path
  
  dfs <- list()
  for (file in files.auc) {
    dfs <- append(dfs, list(fread(file)))
  }
  
  df.per.scenario <- rbindlist(dfs)
  df.per.scenario\$scenario <- scenario
  
  df <- rbind(df, df.per.scenario)
}

df <- df[method %in% c("deseq2", "dream")]

# Introduce new column for plotting to avoid long scenario names
replacement_vector <- c("dataset-ub-cells" = FALSE, "dataset-ub-cells_pb-fixed-effect" = TRUE)
df[, fixed_effect := replacement_vector[scenario]]

p <- ggplot(df, aes(x = method, y = auc, fill = fixed_effect)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("DESeq2", "DREAM")) +
  ylim(0, 1) +
  labs(y = "AUPRC", fill = "Cell count as\nfixed effect") +
  scale_fill_manual(values = c("#FF7F00", "#1F78B4")) +
  theme_cowplot() +
  theme(axis.title.x = element_blank())

ggsave(p, filename = "Fig_S02.png", width = 2650, height = 2000, units = 'px', dpi = 300)
