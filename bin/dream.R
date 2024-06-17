#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(variancePartition)
  library(edgeR)
  library(anndataR)
  library(argparser)
  library(data.table)
  library(SingleCellExperiment)
})

set.seed(0)

p <- arg_parser("DREAM Script")

p <- add_argument(p, "--input", help = "path to input h5ad")
p <- add_argument(p, "--output", help = "path to output tsv")
p <- add_argument(p, "--scenario", help = "simulated scenario")
p <- add_argument(p, "--threads", help = "number of threads for parallel work", default=4)

argv <- parse_args(p)

input <- argv$input
output <- argv$output
scenario <- argv$scenario
threads <- argv$threads

print(input)
print(output)
print(scenario)
print(threads)


data <- read_h5ad(input, to='SingleCellExperiment')

isexpr <- rowSums(edgeR::cpm(SummarizedExperiment::assays(data)$X) > 0.1) >= 5

dge <- DGEList(SummarizedExperiment::assays(data)$X[isexpr,])
dge <- calcNormFactors(dge)

param <- SnowParam(as.numeric(threads), "SOCK", progressbar = T)

colData(data)$log.n_obs <- log(colData(data)$n_obs)

if (scenario %in% c("atlas", "atlas_hvg", "atlas-less-de")) {
  form <- ~ Condition + (1|Batch)
} else if (scenario %in% c("dataset", "dataset_hvg", "dataset-less-de")) {
  form <- ~ Condition
} else if (scenario %in% c("atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de")) {
  form <- ~ Condition + (1|Batch)
} else if (scenario == "atlas-negative") {
  form <- ~ Condition + (1|Batch)
} else if (scenario == "kang2018") {
  form <- ~ Condition + (1|Batch)
} else if (scenario %in% c("dataset-ub-cells", "dataset-ub-cells_hvg", "dataset-ub-cells-less-de")) {
  form <- ~ Condition
} else if (scenario == "dataset-ub-cells_pb-fixed-effect") {
  form <- ~ Condition + log.n_obs
}

vobjDream <- voomWithDreamWeights(dge, form, colData(data))

fitmm <- dream(vobjDream, form, colData(data))
fitmm <- eBayes(fitmm)

dt <- data.table(limma::topTable(fitmm, number=Inf), keep.rownames = "gene")[order(gene)]

fwrite(dt, file = output, sep = "\t")
