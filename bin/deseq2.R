#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(anndataR)
  library(argparser)
})

set.seed(0)

# Parse arguments
p <- arg_parser("DESeq2 Script")

p <- add_argument(p, "--input", help = "path to input h5ad")
p <- add_argument(p, "--output", help = "path to output tsv")
p <- add_argument(p, "--scenario", help = "simulation scenario, can be 'atlas', 'dataset' and 'paired'")

argv <- parse_args(p)

input <- argv$input
output <- argv$output
scenario <- argv$scenario

print(input)
print(output)
print(scenario)

print("Loading data...")

sce <- read_h5ad(input, to='SingleCellExperiment')

print("Data loaded")

if ("Batch" %in% colnames(colData(sce))) {
  colData(sce)$Batch <- as.factor(colData(sce)$Batch)
}
colData(sce)$Sample <- as.factor(colData(sce)$Sample)

print("Converted colData")

counts <- assays(sce)$X
colnames(counts) <- colData(sce)[, 'Sample']
colData <- as.data.frame(colData(sce))
rownames(colData) <- colData[, "Sample"]
colData$log.n_obs <- log(colData$n_obs)

if (scenario %in% c("atlas", "atlas_hvg", "atlas-less-de")) {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ Batch + Condition)
} else if (scenario %in% c("dataset", "dataset_hvg", "dataset-less-de")) {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ Condition)
} else if (scenario %in% c("atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de")) {
   dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ Batch + Condition)
} else if (scenario == "atlas-negative") {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ Batch + Condition)
} else if (scenario == "kang2018") {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ Batch + Condition)
} else if (scenario %in% c("dataset-ub-cells", "dataset-ub-cells_hvg", "dataset-ub-cells-less-de")) {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ Condition)
} else if (scenario == "dataset-ub-cells_pb-fixed-effect") {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ log.n_obs + Condition)
} else if (scenario == "luca") {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ Condition)
}


dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)

write.table(res, file = output, sep = "\t")
