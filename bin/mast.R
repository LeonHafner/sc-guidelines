#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MAST)
  library(anndataR)
  library(argparser)
})

# Parse arguments
p <- arg_parser("MAST Script")

p <- add_argument(p, "--input", help = "path to input h5ad")
p <- add_argument(p, "--output", help = "path to output h5ad")
p <- add_argument(p, "--scenario", help = "simulation scenario, can be 'atlas', 'dataset' and 'paired'")

argv <- parse_args(p)

input <- argv$input
output <- argv$output
scenario <- argv$scenario

print(input)
print(output)
print(scenario)


data <- read_h5ad(input, to='SingleCellExperiment')

names(assays(data))[1] <- "counts"

exprsArray <- as.matrix(counts(data))
cData <- colData(data)
fData <- data.frame(Gene = rownames(exprsArray))

exprsArray.log <- log2(exprsArray + 1)

sca <- FromMatrix(exprsArray = exprsArray.log, cData = cData, fData = fData, check_sanity = F)

colData(sca)$Batch <- factor(colData(sca)$Batch)
colData(sca)$Condition <- factor(colData(sca)$Condition)
colData(sca)$Sample <- factor(colData(sca)$Sample)

if (scenario %in% c("atlas", "atlas_hvg", "atlas-less-de")) {
    zlmCond <- zlm(~ Condition + (1 | Batch) + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario %in% c("dataset", "dataset_hvg", "dataset-less-de")) {
    zlmCond <- zlm(~ Condition + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario %in% c("atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de")) {
    zlmCond <- zlm(~ Condition + Batch + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario == "atlas-negative") {
  zlmCond <- zlm(~ Condition + (1 | Batch) + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario == "kang2018") {
  zlmCond <- zlm(~ Condition + Batch + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario %in% c("dataset-ub-cells", "dataset-ub-cells_hvg", "dataset-ub-cells-less-de")) {
  zlmCond <- zlm(~ Condition + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
}

print("Finished zlm")

summaryCond <- MAST::summary(zlmCond, doLRT = 'ConditionCondition2')
summaryDt <- summaryCond$datatable

write.table(summaryDt, file = output, sep = '\t', row.names = FALSE)
