#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)


prc_path_string <- "${prc_path_string}"
auc_meta_string <- "${auc_meta_string}"
auc_path_string <- "${auc_path_string}"

# Parse PRC input path string into dataframe with columns 'scenario', 'method' and 'path'
parse_prc_path_string_into_dataframe <- function(path_string) {
  paths <- strsplit(path_string, ";")[[1]]
  
  df <- data.frame(scenario = character(0), method = character(0), path = character(0))
  
  for (p in paths) {
    path_parts <- unlist(strsplit(p, "/"))
    filename <- gsub(".tsv", "", path_parts[length(path_parts)])
    filename_parts <- strsplit(filename, "_")[[1]]
    
    scenario <- filename_parts[2]
    method <- filename_parts[4]
    path <- p
    df <- rbind(df, data.frame(scenario, method, path))
  }
  return(df)
}

# Split each individual AUC meta object into 'scenario' and 'run'
auc_parse_meta <- function(auc_meta) {
  parts <- strsplit(auc_meta, ", ")[[1]]
  scenario <- sub("scenario:", "", parts[1])
  run <- sub("run:", "", parts[2])
  return(list("scenario" = scenario, "run" = run))
}

# Convert PRC paths to dataframe
prc_paths_df <- parse_prc_path_string_into_dataframe(prc_path_string)
prc_paths_df <- prc_paths_df[prc_paths_df\$method == "distinct",]

# Convert AUC meta and paths to dataframe
auc_meta_components <- strsplit(gsub("\\\\[|\\\\]", "", auc_meta_string), ";")[[1]]
auc_path_components <- strsplit(auc_path_string, ";")[[1]]

auc_data <- lapply(1:length(auc_meta_components), function(i) {
  auc_meta_parsed <- auc_parse_meta(auc_meta_components[i])
  list("scenario" = auc_meta_parsed\$scenario,
       "run" = auc_meta_parsed\$run,
       "path" = auc_path_components[i])
})
auc_paths_df <- do.call(rbind, lapply(auc_data, as.data.frame))

# Read PRC paths into prc data table
prc <- data.table(precision = numeric(), recall = numeric(), method = character(), scenario = character())
for (i in 1:nrow(prc_paths_df)) {
  prc_per_method <- fread(prc_paths_df[i,]\$path)
  prc_per_method\$method <- prc_paths_df[i,]\$method
  prc_per_method\$scenario <- prc_paths_df[i,]\$scenario
  prc <- rbind(prc, prc_per_method)
}

# Read AUC paths into auc data table
auc <- data.table(method = character(), auc = numeric(), scenario = character())
for (scenario in unique(auc_paths_df\$scenario)) {
  files.auc <- auc_paths_df[auc_paths_df\$scenario == scenario,]\$path
  
  dfs <- list()
  for (file in files.auc) {
    dfs <- append(dfs, list(fread(file)))
  }
  
  df.per.scenario <- rbindlist(dfs)
  df.per.scenario\$scenario <- scenario
  
auc <- rbind(auc, df.per.scenario)
}
auc <- auc[method == "distinct"]

color.code <- data.table(scenario = c("atlas", "atlas-ub-conditions", "dataset", "dataset-ub-cells"), 
                         color = c(1:4), 
                         scenario_legend = c("Atlas", "Atlas (unbalanced conditions)", "Dataset", "Dataset (unbalanced cells)"))

auc.ordered <- auc[,.(median_auc = median(auc)), by = "scenario"][order(median_auc, decreasing = T)]

color.code.auc <- merge(auc.ordered, color.code, by = "scenario", sort = F)

prc\$scenario <- factor(prc\$scenario, levels = color.code.auc\$scenario)
auc\$scenario <- factor(auc\$scenario, levels = color.code.auc\$scenario)


plot.prc <- ggplot(prc, aes(x = recall, y = precision, color = scenario)) +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") +
  geom_line() +
  ylim(0, 1) +
  scale_color_okabe_ito(order = color.code.auc\$color, labels = color.code.auc\$scenario_legend) +
  labs(x = "Recall", y = "Precision", color = "Scenario") +
  theme_cowplot()

plot.auc <- ggplot(auc, aes(x = scenario, y = auc, color = scenario)) +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed") + 
  geom_boxplot() +
  ylim(0, 1) +
  labs(color = "Scenario", y = "AUPRC") +
  scale_color_okabe_ito(order = color.code.auc\$color, labels = color.code.auc\$scenario_legend) +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

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

ggsave(p, filename = "Fig_S07_run${meta_prc.run}.png", width = 3606, height = 1950, units = "px", dpi = 400)





