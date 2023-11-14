suppressPackageStartupMessages({
  library("splatter")
  library("zellkonverter")
})

base.path <- 'path/to/base-dir'

source('splimp.R')

# Erhalte die Liste der Argumente
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
nGenes <- as.integer(args[2])
nCells <- as.integer(args[3])
mode <- args[4]

output <- file.path(base.path, 'data', mode, '1_sims', filename)

vcf <- mockVCF(n.samples = 10)

params <- newSplatPopParams(nGenes = nGenes,
                            #pop.mean.shape = 0.34,
                            #pop.mean.rate = 0.008,
                            eqtl.ES.shape = 0,
                            eqtl.ES.rate = 0,
                            #mean.shape = 0.6,
                            #mean.rate = 0.3,
                            #lib.loc = 11,
                            #lib.scale = 0.2,
                            out.prob = 0,
                            #out.facLoc = 4,
                            #out.facScale = 0.5,
                            #bcv.common = 0.1,
                            #bcv.df = 60,
                            #dropout.mid = 0,
                            #dropout.shape = -1,
                            #similarity.scale = 1,
                            #pop.cv.bins = 10,
                            #pop.quant.norm = TRUE,
                            eqtl.n = 0,
                            #eqtl.dist = 1e6,
                            #eqtl.maf.min = 0.05,
                            #eqtl.maf.max = 0.5,
                            #eqtl.coreg = 0,
                            eqtl.group.specific = 0,
                            eqtl.condition.specific = 0,
                            #nCells.sample = FALSE,
                            #nCells.shape = 1.5,
                            #nCells.rate = 0.015,
                            batch.size = 10,
                            batchCells = c(nCells / 10),
                            #batch.facLoc = c(0.0),
                            #batch.facScale = c(0.0),
                            #group.prob = 1
                            #de.prob = 0.1,
                            #de.downProb = 0.5,
                            #de.facLoc = 0.1,
                            #de.facScale = 0.4,
                            condition.prob = c(0.5, 0.5),
                            cde.prob = 0.025,
                            cde.downProb = 0.5,
                            cde.facLoc = 2.5,
                            cde.facScale = 0.4
                            )

sim <- splimpSimulate(vcf = vcf, params = params)

assays(sim)$BCV <- NULL
assays(sim)$BaseCellMeans <- NULL
assays(sim)$BatchCellMeans <- NULL
assays(sim)$CellMeans <- NULL
assays(sim)$TrueCounts <- NULL

writeH5AD(sim, output, X_name = 'counts')
