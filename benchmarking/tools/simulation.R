suppressPackageStartupMessages({
  library("splatter")
  library("zellkonverter")
  library("argparser")
})


# Parse arguments
p <- arg_parser("Simulation Script")

p <- add_argument(p, "--output", help = "path for output h5ad")
p <- add_argument(p, "--scenario", help = "simulation scenario, can be 'atlas', 'dataset' and 'paired'")

argv <- parse_args(p)

output <- argv$output
scenario <- argv$scenario

# Simulation
source("tools/splimp.R")


if (scenario == "atlas") {
    vcf <- mockVCF(n.samples = 15)
    params <- newSplatPopParams(nGenes = 20000,
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
} else if (scenario == "dataset") {
    vcf <- mockVCF(n.samples = 10)
    params <- newSplatPopParams(nGenes = 20000,
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
} else if (scenario == "atlas-ub-conditions") {
  vcf <- mockVCF(n.samples = 32)
  params <- newSplatPopParams(nGenes = 20000,
                              eqtl.ES.shape = 0,
                              eqtl.ES.rate = 0,
                              out.prob = 0,
                              eqtl.n = 0,
                              eqtl.group.specific = 0,
                              eqtl.condition.specific = 0,
                              batch.size = 10,
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
  params <- newSplatPopParams(nGenes = 20000,
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
} else if (scenario == "dataset-ub-cells") {
  vcf <- mockVCF(n.samples = 10)
  params <- newSplatPopParams(nGenes = 20000,
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
} else if (scenario == "atlas_less_de") {
  vcf <- mockVCF(n.samples = 15)
  params <- newSplatPopParams(nGenes = 20000,
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
} else if (scenario == "atlas-ub-conditions_less_de") {
  vcf <- mockVCF(n.samples = 32)
  params <- newSplatPopParams(nGenes = 20000,
                              eqtl.ES.shape = 0,
                              eqtl.ES.rate = 0,
                              out.prob = 0,
                              eqtl.n = 0,
                              eqtl.group.specific = 0,
                              eqtl.condition.specific = 0,
                              batch.size = 10,
                              batchCells = c(250, 250),
                              batch.facLoc = c(0.05, 0.8),
                              batch.facScale = c(0.1, 0.1),
                              condition.prob = c(0.5, 0.5),
                              cde.prob = 0.0025,
                              cde.downProb = 0.5,
                              cde.facLoc = 2.5,
                              cde.facScale = 0.4
  )
} else if (scenario == "dataset_less_de") {
  vcf <- mockVCF(n.samples = 10)
  params <- newSplatPopParams(nGenes = 20000,
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
} else if (scenario == "dataset-ub-cells_less_de") {
  vcf <- mockVCF(n.samples = 10)
  params <- newSplatPopParams(nGenes = 20000,
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
}

sim <- splimpSimulate(vcf = vcf, params = params)

assays(sim)$BCV <- NULL
assays(sim)$BaseCellMeans <- NULL
assays(sim)$BatchCellMeans <- NULL
assays(sim)$CellMeans <- NULL
assays(sim)$TrueCounts <- NULL

writeH5AD(sim, output, X_name = 'counts')
