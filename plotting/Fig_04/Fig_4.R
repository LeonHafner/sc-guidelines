library(data.table)
library(ggplot2)
library(cowplot)

setwd('/nfs/home/students/l.hafner/sc-guidelines/plotting/Fig_4')

data <- rbindlist(lapply(paste0('correlations_run', 0:9, '.tsv'), FUN = fread))

color.intra <- '#ff7f00'
color.inter <- '#1f78b4'

p1 <- ggplot(data, aes(x = run, y = correlation, fill = type)) +
        geom_boxplot() +
        labs(x = 'Simulation', y = 'Correlation (spearman)', fill = "Type of\ncorrelation") +
        theme(axis.text = element_text(size=10),
              axis.title = element_text(size=12),
              legend.text = element_text(size=10),
              legend.title = element_text(size=12)) +
        scale_x_discrete(labels = paste('Sim', 1:10, sep = '')) +
        scale_fill_manual(values = c(color.inter, color.intra)) +
        theme_cowplot(14)

ggsave(p1, file = "Fig_4.png", width = 2480, height = 1351, units = 'px')
