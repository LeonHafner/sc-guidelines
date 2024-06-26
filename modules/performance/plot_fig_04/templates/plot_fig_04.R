#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)
library(pracma)


prc_path_string <- "${prc_path_string}"

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

prc_paths_df <- parse_prc_path_string_into_dataframe(prc_path_string)

# Read PRC paths into prc data table
prc <- data.table(precision = numeric(), recall = numeric(), method = character(), scenario = character())
for (i in 1:nrow(prc_paths_df)) {
  prc_per_method <- fread(prc_paths_df[i,]\$path)
  prc_per_method\$method <- prc_paths_df[i,]\$method
  prc_per_method\$scenario <- prc_paths_df[i,]\$scenario
  prc <- rbind(prc, prc_per_method)
}

prc <- na.omit(prc)
prc <- prc[(recall != 0 | precision != 0)]


# Unified color coding of the methods
color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation-test", "scvi", "distinct"), 
                         color = c(1:7),
                         hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "distinct"))

# Factorize method column to fit colors
prc\$method <- factor(prc\$method, levels = color.code\$method)


# Compute AUC of the curves
auc.deseq2 <- round(trapz(prc[method == 'deseq2', recall], prc[method == 'deseq2', precision]), 3)
auc.mast <- round(trapz(prc[method == 'mast', recall], prc[method == 'mast', precision]), 3)
auc.hb <- round(trapz(prc[method == 'hierarchical-bootstrapping', recall], prc[method == 'hierarchical-bootstrapping', precision]), 3)
auc.perm <- round(trapz(prc[method == 'permutation-test', recall], prc[method == 'permutation-test', precision]), 3)
auc.scvi <- round(trapz(prc[method == 'scvi', recall], prc[method == 'scvi', precision]), 3)
auc.dream <- round(trapz(prc[method == 'dream', recall], prc[method == 'dream', precision]), 3)
auc.distinct <- round(trapz(prc[method == 'distinct', recall], prc[method == 'distinct', precision]), 3)

p <- ggplot(prc, aes(x = recall, y = precision, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = tail(prc\$precision, n = 1), linetype = "dashed", color = "red") +
  labs(x = "Recall", y = "Precision", color = 'Method') +
  annotate("text", x = 0.75, y = 0.77, label = "AUC:", size = 4, color = "black") +
  annotate("text", x = 0.75, y = 0.70, label = paste(auc.deseq2), size = 5, color = color.code[method == "deseq2", hex]) +
  annotate("text", x = 0.75, y = 0.63, label = paste(auc.dream), size = 5, color = color.code[method == "dream", hex]) +
  annotate("text", x = 0.75, y = 0.56, label = paste(auc.hb), size = 5, color = color.code[method == "hierarchical-bootstrapping", hex]) +
  annotate("text", x = 0.75, y = 0.49, label = paste(auc.mast), size = 5, color = color.code[method == "mast", hex]) +
  annotate("text", x = 0.75, y = 0.42, label = paste(auc.perm), size = 5, color = color.code[method == "permutation-test", hex]) +
  annotate("text", x = 0.75, y = 0.35, label = paste(auc.scvi), size = 5, color = color.code[method == "scvi", hex]) +
  annotate("text", x = 0.75, y = 0.28, label = paste(auc.distinct), size = 5, color = color.code[method == "distinct", hex]) +
  scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_cowplot(16) +
  theme(axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16))

ggsave(p, filename = "Fig_04.png", width = 2650, height = 2000, units = "px", dpi = 300)