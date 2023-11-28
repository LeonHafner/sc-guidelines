library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggokabeito)

setwd('path/to/base-dir')

results.genes.fixed <- fread('results_genes_fixed_mean.tsv')
results.cells.fixed <- fread('results_cells_fixed_mean.tsv')


# Add time for pseudobulking
results.genes.fixed$deseq2 <- results.genes.fixed$deseq2 + results.genes.fixed$pseudobulk
results.genes.fixed$permutation <- results.genes.fixed$permutation + results.genes.fixed$pseudobulk
results.genes.fixed$dream <- results.genes.fixed$dream + results.genes.fixed$pseudobulk

results.cells.fixed$deseq2 <- results.cells.fixed$deseq2 + results.cells.fixed$pseudobulk
results.cells.fixed$permutation <- results.cells.fixed$permutation + results.cells.fixed$pseudobulk
results.cells.fixed$dream <- results.cells.fixed$dream + results.cells.fixed$pseudobulk


# Unified color coding of the methods
color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation", "scvi", "distinct"), 
                         color = c(1:7), 
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "distinct"))


# Melt and filter data for the right columns
results.cells.fixed <- results.cells.fixed %>%
                        melt(id.vars = c('n_genes', 'n_cells'),
                             value.name = "time",
                             variable.name = "process") %>%
                        filter(!process %in% c('de-genes', 'pseudobulk', 'sim', 'preprocess'))

results.genes.fixed <- results.genes.fixed %>%
                        melt(id.vars = c('n_genes', 'n_cells'),
                             value.name = "time",
                             variable.name = "process") %>%
                        filter(!process %in% c('de-genes', 'pseudobulk', 'sim', 'preprocess'))


# Make factor to change order of legend
results.cells.fixed$process <- factor(results.cells.fixed$process, levels = color.code$method)
results.genes.fixed$process <- factor(results.genes.fixed$process, levels = color.code$method)
  
  
p1 <-   ggplot(results.cells.fixed, aes(x = n_genes, y = time, color = process)) +
          geom_line(linewidth=0.8) +
          labs(x = '#Genes', y = 'Time (s)', color = 'Method') +
          scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
          theme_cowplot(14) +
          background_grid(major = "xy",
                          minor = "none",
                          size.major = 0.2,
                          size.minor = 0.2)


p2 <- ggplot(results.genes.fixed, aes(x = n_cells, y = time, color = process)) +
        geom_line(linewidth=0.8) +
        labs(x = '#Cells', y = 'Time (s)', color = 'Method') +
        scale_color_okabe_ito(order = color.code$color, labels = color.code$method_legend) +
        theme_cowplot(14) +
        background_grid(major = "xy",
                        minor = "none",
                        size.major = 0.2,
                        size.minor = 0.2)


legend.double.plot <- get_legend(p1 + theme(legend.margin = margin(0, 0, 20, 20)))

double.plot <- plot_grid(p1 + theme(legend.position = "none"),
                         p2 + theme(legend.position = "none"),
                         labels = "AUTO",
                         label_size = 20,
                         ncol = 1,
                         label_x = 0, label_y = 1,
                         hjust = -0.1, vjust = 1.1,
                         align = "v")


runtime_benchmark <- plot_grid(double.plot, legend.double.plot, rel_widths = c(2.5, 1))


ggsave(runtime_benchmark, width = 2370, height = 3194, units = "px", dpi = 300, filename = "Fig_6.png")
