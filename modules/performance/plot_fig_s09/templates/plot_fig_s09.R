#!/usr/bin/env Rscript

library(data.table)
library(anndataR)
library(SingleCellExperiment)
library(dplyr)
library(cowplot)
library(ggplot2)

sce <- read_h5ad("${input_anndata}", to = "SingleCellExperiment")
run <- "${meta.run}"

colData <- colData(sce) %>%
  as.data.table() %>%
  mutate(Sample = as.factor(rep(seq(0, 19))), Condition = factor(gsub("Condition", "", Condition), levels = c(1 , 2)), n_obs = NULL)

rowData <- as.data.table(rowData(sce), keep.rownames = "gene")
rowData[, foldchange := pmax(ConditionDE.Condition1 / ConditionDE.Condition2, ConditionDE.Condition2 / ConditionDE.Condition1)]
top_rowData <- rowData[order(foldchange, decreasing = T), .(gene, foldchange)][foldchange > 1,][1:8]
low_rowData <- rowData[foldchange == 1, .(gene, foldchange)][sample(.N, 8, replace = FALSE)]

data <- cbind(colData, t(assays(sce)\$X[c(top_rowData\$gene, low_rowData\$gene),]))

data <- data %>%
  melt(id.vars = c("Batch", "Sample", "Condition"), variable.name = "gene", value.name = "count")

data <- merge(data, rbindlist(list(top_rowData, low_rowData)), by = "gene")
data <- mutate(data, gene_fc = paste0(gene, " ", "[FC=", round(foldchange, 2), "]"))

data <- data[order(foldchange, gene_fc)]

data\$gene_fc <- factor(data\$gene_fc, levels = unique(data\$gene_fc))

p <- ggplot(data, aes(x = Batch, y = count, fill=Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#1F78B4", "#FF7F00")) +
  labs(y = "Count") +
  scale_y_log10() +
  facet_wrap(~ gene_fc, scales="free") +
  theme_cowplot() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 8))

ggsave(p, filename = "Fig_S09_run${meta.run}.png", width = 2480, height = 3100, units = "px", dpi=300)
