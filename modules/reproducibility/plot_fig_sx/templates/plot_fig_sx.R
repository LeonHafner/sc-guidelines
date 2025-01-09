#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)

data_dir <- "."
N_fixed <- 100

# 1) Read all files
files <- list.files(data_dir, pattern = ".tsv\$", full.names = TRUE)
dt_list <- lapply(files, fread)
names(dt_list) <- basename(files)

# 2) Get the intersection of all genes across data.tables
all_genes <- Reduce(intersect, lapply(dt_list, function(dt) dt\$V1))
all_genes <- sort(all_genes)

# 3) Subset each data.table to keep only the intersecting genes
dt_list <- lapply(dt_list, function(dt) dt[dt\$V1 %in% all_genes])

# 4) Identify methods (assumes first column is gene name, subsequent columns are methods)
methods <- colnames(dt_list[[1]])[-1]
num_datasets <- length(dt_list)
num_genes <- length(all_genes)

compute_jaccards_for_N <- function(N, methods, dt_list) {
  method_jaccards <- list()
  
  for (method in methods) {
    # Collect top-N sets for this method across all data.tables
    top_sets <- list()
    for (dt in dt_list) {
      dt_sorted <- dt[order(get(method))]
      top_genes <- dt_sorted\$V1[1:N]
      top_sets[[length(top_sets) + 1]] <- top_genes
    }
    
    # Compute pairwise Jaccard among all datasets
    combs <- combn(seq_along(top_sets), 2, simplify = FALSE)
    jaccards <- numeric(length(combs))
    for (i in seq_along(combs)) {
      s1 <- top_sets[[combs[[i]][1]]]
      s2 <- top_sets[[combs[[i]][2]]]
      intersection_count <- length(intersect(s1, s2))
      union_count <- length(union(s1, s2))
      jaccards[i] <- if (union_count > 0) intersection_count / union_count else 0
    }
    
    method_jaccards[[method]] <- jaccards
  }
  
  return(method_jaccards)
}

fixed_jaccards <- compute_jaccards_for_N(N_fixed, methods, dt_list)

boxplot_data <- data.frame(Method = character(), Jaccard = numeric())
for (method in names(fixed_jaccards)) {
  jaccs <- fixed_jaccards[[method]]
  boxplot_data <- rbind(boxplot_data, 
                        data.frame(Method = method, Jaccard = jaccs))
}



color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation-test", "scvi", "distinct", "ttest"), 
                         color = c(1, 2, 3, 4, 5, 6, 7, 9),
                         hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"),
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "distinct", "t-test"))

boxplot_data\$Method <- factor(boxplot_data\$Method, levels = color.code\$method)

p <- ggplot(boxplot_data, aes(x = Method, y = Jaccard, fill=Method)) +
  geom_boxplot() +
  scale_fill_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
  scale_x_discrete(labels = color.code\$method_legend) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p, filename = "Fig_SX.png", width = 2048, height = 1500, units = "px")
