#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ggokabeito)
  library(data.table)
  library(pracma)
})

path_string <- "${path_string}"


parse_path_string_into_dataframe <- function(path_string) {
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

plot_prc_of_scenario <- function(df, scenario) {
  files_for_scenario <- df[df\$scenario == scenario,]
  prc <- data.table(precision = numeric(), recall = numeric(), method = character())
  
  for (i in 1:nrow(files_for_scenario)) {
    prc_per_method <- fread(files_for_scenario[i,]\$path)
    prc_per_method\$method <- files_for_scenario[i,]\$method
    prc <- rbind(prc, prc_per_method)
  }
  
  prc\$method <- factor(prc\$method, levels = c("deseq2", "hierarchical-bootstrapping", "mast", "permutation-test", "distinct", "scvi", "dream", "ttest"))
  
  prc <- na.omit(prc)
  prc <- prc[(recall != 0 | precision != 0)]
  
  color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation-test", "scvi", "ttest"), 
                           color = c(1:7), 
                           hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                           method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "t-test"))
  
  # Remove distinct as those are bad values
  prc <- prc[method != "distinct"]
  
  
  # Factorize method column to fit colors
  prc\$method <- factor(prc\$method, levels = color.code\$method)
  
  # Compute AUC of the curves
  auc.deseq2 <- round(trapz(prc[method == 'deseq2', recall], prc[method == 'deseq2', precision]), 3)
  auc.mast <- round(trapz(prc[method == 'mast', recall], prc[method == 'mast', precision]), 3)
  auc.hb <- round(trapz(prc[method == 'hierarchical-bootstrapping', recall], prc[method == 'hierarchical-bootstrapping', precision]), 3)
  auc.perm <- round(trapz(prc[method == 'permutation-test', recall], prc[method == 'permutation-test', precision]), 3)
  auc.scvi <- round(trapz(prc[method == 'scvi', recall], prc[method == 'scvi', precision]), 3)
  auc.dream <- round(trapz(prc[method == 'dream', recall], prc[method == 'dream', precision]), 3)
  auc.ttest <- round(trapz(prc[method == 'ttest', recall], prc[method == 'ttest', precision]), 3)
  
  
  p <- ggplot(prc, aes(x = recall, y = precision, color = method)) +
    geom_line(linewidth = 0.8) +
    xlim(0, 1) +
    ylim(0, 1) +
    geom_hline(yintercept = tail(prc\$precision, n = 1), linetype = "dashed", color = "red") +
    labs(x = "Recall", y = "Precision", color = 'Method') +
    annotate("text", x = 0.1, y = 0.63, label = "AUC:", size = 4, color = "black") +
    annotate("text", x = 0.1, y = 0.56, label = paste(auc.deseq2), size = 4, color = color.code[method == "deseq2", hex]) +
    annotate("text", x = 0.1, y = 0.49, label = paste(auc.dream), size = 4, color = color.code[method == "dream", hex]) +
    annotate("text", x = 0.1, y = 0.42, label = paste(auc.hb), size = 4, color = color.code[method == "hierarchical-bootstrapping", hex]) +
    annotate("text", x = 0.1, y = 0.35, label = paste(auc.mast), size = 4, color = color.code[method == "mast", hex]) +
    annotate("text", x = 0.1, y = 0.28, label = paste(auc.perm), size = 4, color = color.code[method == "permutation-test", hex]) +
    annotate("text", x = 0.1, y = 0.21, label = paste(auc.scvi), size = 4, color = color.code[method == "scvi", hex]) +
    annotate("text", x = 0.1, y = 0.14, label = paste(auc.ttest), size = 4, color = color.code[method == "ttest", hex]) +
    scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
    theme_cowplot() +
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 10))

  return(p)
}

df <- parse_path_string_into_dataframe(path_string)

p.dataset <- plot_prc_of_scenario(df, "dataset")
p.dataset.ub <- plot_prc_of_scenario(df, "dataset-ub-cells")
p.atlas <- plot_prc_of_scenario(df, "atlas")
p.atlas.ub <- plot_prc_of_scenario(df, "atlas-ub-conditions")

p <- plot_grid(
  plot_grid(p.atlas + theme(axis.title.x = element_blank(),
                            legend.position = "none"),
            p.atlas.ub + theme(axis.title.y = element_blank(),
                               axis.title.x = element_blank(),
                               legend.position = "none"),
            p.dataset + theme(legend.position = "none"),
            p.dataset.ub + theme(axis.title.y = element_blank(),
                                 legend.position = "none"),
            labels = "AUTO",
            label_size = 20,
            ncol = 2,
            label_x = 0,
            label_y = 1,
            hjust = -0.1,
            vjust = 1.1,
            align = "hv"),
  get_legend(p.atlas +
               theme(legend.text = element_text(size = 11),
                     legend.title = element_text(size = 15),
                     legend.spacing.x = unit(0.5, "cm"),
                     legend.title.align = 0.5,
                     legend.position=c(0.3, 0.6)) +
               guides(color = guide_legend(ncol = 3))),
  nrow = 2,
  rel_heights = c(1, 0.3),
  vjust = -10
)

ggsave(p, width = 3606, height = 4096, units = "px", dpi = 400, filename = "Fig_S07_run${meta.run}.png")











