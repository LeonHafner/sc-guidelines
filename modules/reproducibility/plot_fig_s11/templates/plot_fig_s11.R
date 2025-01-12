#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)
library(gtools)

data_dir <- "."
files <- list.files(data_dir, pattern = ".tsv\$", full.names = TRUE)
files <- mixedsort(files)

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

results <- list()

for (N in seq_len(num_genes)) {
  method_jaccards <- compute_jaccards_for_N(N, methods, dt_list)
  
  for (method in names(method_jaccards)) {
    jacc_vector <- method_jaccards[[method]]
    avg_jacc <- if (length(jacc_vector) > 0) mean(jacc_vector) else 0
    results[[length(results) + 1]] <- data.frame(N = N, Method = method,Avg_Jaccard = avg_jacc)
  }
  if (N %% 10 == 0) {
    print(paste(N, "/", num_genes))
  }
}

line_df <- do.call(rbind, results)

color.code <- data.table(method = c("deseq2", "dream", "hierarchical-bootstrapping", "mast", "permutation-test", "scvi", "distinct", "scdd", "ttest"), 
                         color = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                         hex = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000"),
                         method_legend = c("DESeq2", "DREAM", "Hierarchical\nBootstrapping", "MAST", "Permutation\nTest", "scVI", "distinct", "scDD", "t-test"))

line_df\$Method <- factor(line_df\$Method, levels = color.code\$method)

p <- ggplot(line_df, aes(x = N, y = Avg_Jaccard, color = Method)) +
  geom_line() +
  scale_color_okabe_ito(order = color.code\$color, labels = color.code\$method_legend) +
  labs(x = "First N genes considered as DE", y = "Avg. Jaccard Index") +
  theme(legend.position = "right") +
  theme_cowplot()

ggsave(p, filename = "Fig_S11.png", width = 2048, height = 1800, units = "px")  
