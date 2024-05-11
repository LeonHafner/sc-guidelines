#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(distinct)
  library(anndataR)
  library(scuttle)
  library(argparser)
})

# Parse arguments
p <- arg_parser("Distinct Script")

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

#Rename assay
names(assays(data))[1] <- "counts"

data.norm <- logNormCounts(data)

colData(data.norm)$Group <- "Group1"
colData <- as.data.frame(colData(data.norm))
unique.samples <- unique(colData[,c("Batch", "Sample", "Condition")])

batch <- unique.samples$Batch
sample <- unique.samples$Sample
condition <- unique.samples$Condition

if (scenario %in% c("atlas", "atlas_hvg", "atlas-less-de")) {
    design <- model.matrix(~ condition + batch)
} else if (scenario %in% c("dataset", "dataset_hvg", "dataset-less-de")) {
    design <- model.matrix(~ condition)
} else if (scenario %in% c("atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de")) {
   design <- model.matrix(~ condition + batch)
} else if (scenario == "atlas-negative") {
  design <- model.matrix(~ condition + batch)
} else if (scenario == "kang2018"){
  design <- model.matrix(~ condition + batch)
} else if (scenario %in% c("dataset-ub-cells", "dataset-ub-cells_hvg", "dataset-ub-cells-less-de")) {
  design <- model.matrix(~ condition)
}

rownames(design) <- sample

res <- distinct_test(x = data.norm,
                     name_assays_expression = 'logcounts',
                     name_cluster = 'Group',
                     name_sample = 'Sample',
                     design = design,
                     column_to_test = 2,
                     min_non_zero_cells = 20,
                     n_cores = 8)


write.table(res, file = output, sep = '\t', row.names = FALSE)
