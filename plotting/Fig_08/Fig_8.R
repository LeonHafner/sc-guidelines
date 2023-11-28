suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(pracma)
  library(cowplot)
  library(tidyr)
  library(argparser)
})
  
setwd('path/to/base-dir')

methods <- c('mast', 'distinct', 'deseq2', 'permutation', 'hierarchical-bootstrapping', 'scvi', 'dream')

files <- paste0(methods, '.tsv')

prc <- data.table(precision = numeric(), recall = numeric(), method = character())
for (file in files) {
  method <- fread(file)
  method$method <- strsplit(strsplit(file, split = '.', fixed = T)[[1]][1], split = '_', fixed = T)[[1]]
  prc <- rbind(prc, method)
}

prc <- na.omit(prc)
prc <- prc[(recall != 0 | precision != 0)]


# Color coding of the methods
color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation", "scvi", "distinct"), 
                         color = c(1:7), 
                         hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "distinct"))

prc$method <- factor(prc$method, levels = color.code$method)


# Compute AUPRC
auc.deseq2 <- round(trapz(prc[method == 'deseq2', recall], prc[method == 'deseq2', precision]), 3)
auc.mast <- round(trapz(prc[method == 'mast', recall], prc[method == 'mast', precision]), 3)
auc.hb <- round(trapz(prc[method == 'hierarchical-bootstrapping', recall], prc[method == 'hierarchical-bootstrapping', precision]), 3)
auc.perm <- round(trapz(prc[method == 'permutation', recall], prc[method == 'permutation', precision]), 3)
auc.distinct <- round(trapz(prc[method == 'distinct', recall], prc[method == 'distinct', precision]), 3)
auc.scvi <- round(trapz(prc[method == 'scvi', recall], prc[method == 'scvi', precision]), 3)
auc.dream <- round(trapz(prc[method == 'dream', recall], prc[method == 'dream', precision]), 3)


p <- ggplot(prc, aes(x = recall, y = precision, color = method)) +
  geom_line(linewidth = 0.8) +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_hline(yintercept = tail(prc$precision, n = 1), linetype = "dashed", color = "red") +
  labs(x = "Recall", y = "Precision", color = 'Method') +
  annotate("text", x = 0.8, y = 0.93, label = "AUC:", size = 5, color = "black") +
  annotate("text", x = 0.8, y = 0.86, label = paste(auc.deseq2), size = 5, color = color.code[method == "deseq2", hex]) +
  annotate("text", x = 0.8, y = 0.79, label = paste(auc.dream), size = 5, color = color.code[method == "dream", hex]) +
  annotate("text", x = 0.8, y = 0.72, label = paste(auc.hb), size = 5, color = color.code[method == "hierarchical-bootstrapping", hex]) +
  annotate("text", x = 0.8, y = 0.65, label = paste(auc.mast), size = 5, color = color.code[method == "mast", hex]) +
  annotate("text", x = 0.8, y = 0.58, label = paste(auc.perm), size = 5, color = color.code[method == "permutation", hex]) +
  annotate("text", x = 0.8, y = 0.51, label = paste(auc.scvi), size = 5, color = color.code[method == "scvi", hex]) +
  annotate("text", x = 0.8, y = 0.44, label = paste(auc.distinct), size = 5, color = color.code[method == "distinct", hex]) +
  scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
  theme_cowplot(14)



ggsave(p, width = 2650, height = 2000, units = 'px', dpi = 300, filename = "Fig_8.png")

