#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparser)
    library(anndataR)
    library(SingleCellExperiment)
    library(scDD)
})

parser <- arg_parser(description = "Run scDD for single-cell differential expression analysis.")
parser <- add_argument(parser, "--input", type = "character", help = "Input h5ad file (AnnData object).")
parser <- add_argument(parser, "--output", type = "character", help = "Output CSV file for results.")
parser <- add_argument(parser, "--n_cores", type = "integer", default = 1, help = "Number of cores for parallel processing.")

# Parse arguments
args <- parse_args(parser)

print("Test")
sce <- read_h5ad(args$input, to="SingleCellExperiment")
print("Test2")
# Normalization
raw_counts <- assay(sce)
assays(sce)$normcounts <- t(t(raw_counts) / colSums(raw_counts) * 1e6)


prior_param <- list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)

if (inherits(normcounts(sce), "dgCMatrix")) {
    normcounts(sce) <- as.matrix(normcounts(sce))
}

test <- scDD(sce, prior_param = prior_param, testZeroes = T, categorize = F, condition = "Condition")

test.results <- results(test)

output <- test.results[, c("gene", "combined.pvalue")]
colnames(output) <- c("genes", "pvalue")

stopifnot(nrow(output) == nrow(sce))

write.table(output, file = args$output, sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
