suppressPackageStartupMessages({
  library(variancePartition)
  library(edgeR)
  library(BiocParallel)
  library(zellkonverter)
  library(argparser)
  library(data.table)
  library(SingleCellExperiment)
})

p <- arg_parser("DREAM Script")

p <- add_argument(p, "--input", help = "path to input h5ad")
p <- add_argument(p, "--output", help = "path to output tsv")
p <- add_argument(p, "--scenario", help = "simulated scenario")
p <- add_argument(p, "--threads", help = "number of threads for parallel work")

argv <- parse_args(p)

input <- argv$input
output <- argv$output
scenario <- argv$scenario
threads <- argv$threads

print(input)
print(output)
print(scenario)
print(threads)


data <- readH5AD(input)

if (scenario == "kang2018") {
  colnames(colData(data)) <- c("Batch", "Sample", "Condition", "n_obs")
}

isexpr <- rowSums(edgeR::cpm(SummarizedExperiment::assays(data)$X) > 0.1) >= 5

dge <- DGEList(SummarizedExperiment::assays(data)$X[isexpr,])
dge <- calcNormFactors(dge)

param <- SnowParam(as.numeric(threads), "SOCK", progressbar = T)

colData(data)$log.n_obs <- log(colData(data)$n_obs)

if (scenario == "atlas") {
  form <- ~ Condition + (1|Batch)
} else if (scenario == "dataset") {
  form <- ~ Condition
} else if (scenario == "atlas-ub-conditions") {
  form <- ~ Condition + (1|Batch)
} else if (scenario == "atlas-negative") {
  form <- ~ Condition + (1|Batch)
} else if (scenario == "kang2018") {
  form <- ~ Condition + (1|Batch)
} else if (scenario == "dataset-ub-cells") {
  form <- ~ Condition
} else if (scenario == "dataset-ub-cells_pb-fixed-effect") {
  form <- ~ Condition + log.n_obs
}

vobjDream <- voomWithDreamWeights(dge, form, colData(data))

fitmm <- dream(vobjDream, form, colData(data))
fitmm <- eBayes(fitmm)

dt <- data.table(limma::topTable(fitmm, number=Inf), keep.rownames = "gene")[order(gene)]

fwrite(dt, file = output, sep = "\t")
