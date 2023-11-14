suppressPackageStartupMessages({
  library(DESeq2)
  library(zellkonverter)
  library(argparser)
})

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

sce <- readH5AD(input)

print("Data loaded")

if (scenario == "kang2018") {
  colnames(colData(sce)) <- c("Batch", "Sample", "Condition", "n_obs")
}

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

if (scenario == "atlas") {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ Batch + Condition)
} else if (scenario == "dataset") {
    dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ Condition)
} else if (scenario == "atlas-ub-conditions") {
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
} else if (scenario == "dataset-ub-cells") {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ Condition)
} else if (scenario == "dataset-ub-cells_pb-fixed-effect") {
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = colData,
                                design = ~ log.n_obs + Condition)
}


dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)

write.table(res, file = output, sep = "\t")
