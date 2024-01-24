#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library("splatter")
  library("zellkonverter")
  library("argparser")
})


# Parse arguments
p <- arg_parser("Simulation Script")

p <- add_argument(p, "--output", help = "path for output h5ad")
p <- add_argument(p, "--scenario", help = "simulation scenario, can be 'atlas', 'dataset' and 'paired'")
p <- add_argument(p, "--n_genes", help = "number of genes to simulate", type = "integer", default = 5000)

argv <- parse_args(p)

output <- argv$output
scenario <- argv$scenario
n_genes <- argv$n_genes

# Simulation
source("/nfs/data/sc-guidelines/nextflow_test/bin/splimp.R")

getConditionCounts <- function(sim) {
  df <- unique(data.frame(colData(sim)[c("Batch", "Condition", "Sample")]))
  conditions <- expand.grid(Batch = c("Batch1", "Batch2"), 
                            Condition = c("Condition1", "Condition2"))
  apply(conditions, 1, function(cond) {
    nrow(df[df$Batch == cond[1] & df$Condition == cond[2],])
  })
}

assignParams <- function(scenario) {
  if (scenario %in% c("atlas", "atlas_hvg")) {
    vcf <- mockVCF(n.samples = 15)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                batch.size = 5,
                                batchCells = c(250, 250, 250),
                                batch.facLoc = c(0.05, 0.8, 1.45),
                                batch.facScale = c(0.1, 0.1, 0.05),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.025,
                                cde.downProb = 0.5,
                                cde.facLoc = 2.5,
                                cde.facScale = 0.4
    )
  } else if (scenario %in% c("dataset", "dataset_hvg")) {
    vcf <- mockVCF(n.samples = 10)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                batch.size = 10,
                                batchCells = c(250),
                                batch.facLoc = c(0),
                                batch.facScale = c(0),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.025,
                                cde.downProb = 0.5,
                                cde.facLoc = 2.5,
                                cde.facScale = 0.4
    )
  } else if (scenario %in% c("atlas-ub-conditions", "atlas-ub-conditions_hvg")) {
    vcf <- mockVCF(n.samples = 36)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                batch.size = 18,
                                batchCells = c(250, 250),
                                batch.facLoc = c(0.05, 0.8),
                                batch.facScale = c(0.1, 0.1),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.025,
                                cde.downProb = 0.5,
                                cde.facLoc = 2.5,
                                cde.facScale = 0.4
    )
  } else if (scenario == "atlas-negative") {
    vcf <- mockVCF(n.samples = 15)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                batch.size = 5,
                                batchCells = c(250, 250, 250),
                                batch.facLoc = c(0.05, 0.8, 1.45),
                                batch.facScale = c(0.1, 0.1, 0.05),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.025,
                                cde.downProb = 0.5,
                                cde.facLoc = 0,
                                cde.facScale = 0
    )
  } else if (scenario %in% c("dataset-ub-cells", "dataset-ub-cells_hvg")) {
    vcf <- mockVCF(n.samples = 10)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                nCells.sample = TRUE,
                                nCells.shape = 0.8,
                                nCells.rate = 0.0035,
                                batch.size = 10,
                                batchCells = c(250),
                                batch.facLoc = c(0),
                                batch.facScale = c(0),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.025,
                                cde.downProb = 0.5,
                                cde.facLoc = 2.5,
                                cde.facScale = 0.4
    )
  } else if (scenario == "atlas-less-de") {
    vcf <- mockVCF(n.samples = 15)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                batch.size = 5,
                                batchCells = c(250, 250, 250),
                                batch.facLoc = c(0.05, 0.8, 1.45),
                                batch.facScale = c(0.1, 0.1, 0.05),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.0025,
                                cde.downProb = 0.5,
                                cde.facLoc = 2.5,
                                cde.facScale = 0.4
    )
  } else if (scenario == "atlas-ub-conditions-less-de") {
    vcf <- mockVCF(n.samples = 36)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                batch.size = 18,
                                batchCells = c(250, 250),
                                batch.facLoc = c(0.05, 0.8),
                                batch.facScale = c(0.1, 0.1),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.0025,
                                cde.downProb = 0.5,
                                cde.facLoc = 2.5,
                                cde.facScale = 0.4
    )
  } else if (scenario == "dataset-less-de") {
    vcf <- mockVCF(n.samples = 10)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                batch.size = 10,
                                batchCells = c(250),
                                batch.facLoc = c(0),
                                batch.facScale = c(0),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.0025,
                                cde.downProb = 0.5,
                                cde.facLoc = 2.5,
                                cde.facScale = 0.4
    )
  } else if (scenario == "dataset-ub-cells-less-de") {
    vcf <- mockVCF(n.samples = 10)
    params <- newSplatPopParams(nGenes = n_genes,
                                eqtl.ES.shape = 0,
                                eqtl.ES.rate = 0,
                                out.prob = 0,
                                eqtl.n = 0,
                                eqtl.group.specific = 0,
                                eqtl.condition.specific = 0,
                                nCells.sample = TRUE,
                                nCells.shape = 0.8,
                                nCells.rate = 0.0035,
                                batch.size = 10,
                                batchCells = c(250),
                                batch.facLoc = c(0),
                                batch.facScale = c(0),
                                condition.prob = c(0.5, 0.5),
                                cde.prob = 0.0025,
                                cde.downProb = 0.5,
                                cde.facLoc = 2.5,
                                cde.facScale = 0.4
    )
  } else {
    stop("Unknown scenario")
  }
  return(list(params=params, vcf=vcf))
}

assignedParams <- assignParams(scenario)
sim <- splimpSimulate(vcf = assignedParams$vcf, params = assignedParams$params)

if (scenario %in% c("atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de")) {

  conditionCounts <- getConditionCounts(sim)
  
  while (any(conditionCounts < c(9, 1, 1, 9))) {
    assignedParams <- assignParams(scenario)
    sim <- splimpSimulate(vcf = assignedParams$vcf, params = assignedParams$params)

    conditionCounts <- getConditionCounts(sim)
  }
}

writeH5AD(sim, output, X_name = 'counts')