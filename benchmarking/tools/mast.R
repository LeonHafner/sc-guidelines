suppressPackageStartupMessages({
  library(MAST)
  library(scater)
  library(zellkonverter)
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


data <- readH5AD(input)

# Assay is named differently here as this is not simulated
if (scenario == "kang2018") {
  names(assays(data))[1] <- "counts"
  colnames(colData(data)) <- c("Batch", "Condition", "cluster", "cell", "multiplets", "Sample")
}

exprsArray <- as.matrix(counts(data))
cData <- colData(data)
fData <- data.frame(Gene = rownames(exprsArray))

exprsArray.log <- log2(exprsArray + 1)

sca <- FromMatrix(exprsArray = exprsArray.log, cData = cData, fData = fData, check_sanity = F)

colData(sca)$Batch <- factor(colData(sca)$Batch)
colData(sca)$Condition <- factor(colData(sca)$Condition)
colData(sca)$Sample <- factor(colData(sca)$Sample)

if (scenario == "atlas") {
    zlmCond <- zlm(~ Condition + (1 | Batch) + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario == "dataset") {
    zlmCond <- zlm(~ Condition + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario == "atlas-ub-conditions") {
    zlmCond <- zlm(~ Condition + Batch + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario == "atlas-negative") {
  zlmCond <- zlm(~ Condition + (1 | Batch) + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario == "kang2018") {
  zlmCond <- zlm(~ Condition + Batch + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
} else if (scenario == "dataset-ub-cells") {
  zlmCond <- zlm(~ Condition + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
}

if (scenario != "kang2018") {
  summaryCond <- MAST::summary(zlmCond, doLRT = 'ConditionCondition2')
} else {
  summaryCond <- MAST::summary(zlmCond, doLRT = 'Conditionstim')
}
summaryDt <- summaryCond$datatable

write.table(summaryDt, file = output, sep = '\t', row.names = FALSE)
