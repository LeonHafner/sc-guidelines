library(data.table)
library(ggplot2)
library(cowplot)
library(ggokabeito)

path.reg.pb <- "/nfs/home/students/l.hafner/fopra/performance/dataset-ub-cells/data/13_precision-recall"
path.fixed.pb <- "/nfs/home/students/l.hafner/fopra/performance_pb_fixed_effect/dataset-ub-cells/data/13_precision-recall"


df.fe <- data.table(method = character(), auc = numeric(), cell_count_fixed_effect = character())

files.auc.fe <- list.files(path.fixed.pb)
files.auc.fe <- files.auc.fe[startsWith(files.auc.fe, "auc")]

dfs <- list()
for (file in files.auc.fe) {
  dfs <- append(dfs, list(fread(file.path(path.fixed.pb, file))))
}

df.fe <- rbindlist(dfs)
df.fe$cell_count_fixed_effect <- TRUE

df.reg <- data.table(method = character(), auc = numeric(), cell_count_fixed_effect = character())

files.auc.reg <- list.files(path.reg.pb)
files.auc.reg <- files.auc.reg[startsWith(files.auc.reg, "auc")]

dfs <- list()
for (file in files.auc.reg) {
  dfs <- append(dfs, list(fread(file.path(path.reg.pb, file))))
}

df.reg <- rbindlist(dfs)
df.reg$cell_count_fixed_effect <- FALSE

df <- rbind(df.fe, df.reg)

df <- df[method %in% c("dream", "deseq2")]

p <- ggplot(df, aes(x = method, y = auc, fill = cell_count_fixed_effect)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  scale_x_discrete(labels = c("DESeq2", "DREAM")) +
  ylim(0, 1) +
  labs(y = "AUPRC", fill = "Cell count as\nfixed effect") +
  scale_fill_manual(values = c("#FF7F00", "#1F78B4")) +
  theme_cowplot() +
  theme(axis.title.x = element_blank())


ggsave(p, width = 2650, height = 2000, units = 'px', dpi = 300, filename = "/nfs/home/students/l.hafner/fopra/performance_pb_fixed_effect/dataset-ub-cells/data/14_plots/plot.png")
















