#!/usr/bin/env Rscript

# Code is originally from the splatter package
# We extracted it and run it externally to fix some bug that assigned the conditions falsely

suppressPackageStartupMessages({
  library(splatter)
  library(anndataR)
  library(scater)
  library(MAST)
})

splimpSimulate <- function(vcf = mockVCF(n.samples = 20), params = newSplatPopParams(), balancedBatches = TRUE) {
  batch.size <- getParam(params, 'batch.size')
  batchCells <- getParam(params, 'batchCells')
  batch.facLoc <- getParam(params, 'batch.facLoc')
  batch.facScale <- getParam(params, 'batch.facScale')
  message("With Batches:")
  for (param in c('nCells', 'nGenes', 'nBatches', 'batch.size', 'batchCells', 'batch.facLoc', 'batch.facScale')) {
    print(paste0(param, ": ", getParam(params, param)))
  }
  
  params.intern <- params

  params.intern <- setParam(params.intern, "batch.size", 1)
  params.intern <- setParam(params.intern, "batch.facLoc", 0)
  params.intern <- setParam(params.intern, "batch.facScale", 0)
  params.intern <- setParam(params.intern, "batchCells", mean(batchCells))

  message("\nWithout Batches:")
  for (param in c('nCells', 'nGenes', 'nBatches', 'batch.size', 'batchCells', 'batch.facLoc', 'batch.facScale')) {
    print(paste0(param, ": ", getParam(params.intern, param)))
  }
  
  sim <- splatPopSimulate(vcf = vcf, params = params.intern, verbose = T)

  sim <- splimpCleanSCE(sim)
  sim <- splimpDesignBatches(sim, params, balancedBatches)
  sim <- splimpSimBatchEffects(sim, params)
  sim <- splimpSimBatchCellMeans(sim, params)
  sim <- splimpSingleCellMeans(sim)
  sim <- splimpSimBCVMeans(sim, params)
  sim <- splimpSimTrueCounts(sim)
  
  assays(sim)$BCV <- NULL
  assays(sim)$BaseCellMeans <- NULL
  assays(sim)$BatchCellMeans <- NULL
  assays(sim)$CellMeans <- NULL
  assays(sim)$TrueCounts <- NULL

  metadata(sim)$Params <- NULL
  
  return(sim)

}


# Removes unwanted data and data that we are going to simulate from scratch
splimpCleanSCE <- function(sim) {
  assays(sim)$counts <- NULL
  assays(sim)$TrueCounts <- NULL
  assays(sim)$CellMeans <- NULL
  assays(sim)$BCV <- NULL
  assays(sim)$BaseCellMeans <- NULL
  assays(sim)$BatchCellMeans <- NULL
  
  rowData(sim)[, c('chromosome', 'geneStart', 'geneEnd', 'geneMiddle',
                   'eQTL.group', 'eQTL.condition', 'eSNP.ID', 'eSNP.chromosome',
                   'eSNP.loc', 'eSNP.MAF', 'eQTL.EffectSize')] <- NULL
  
  colData(sim)[, c('Batch', 'Group')] <- NULL

  return(sim)
}


# Assign samples to batches
# Pass global params with batch effects here (not params.intern)
# TODO: Change code s.t. there is at least one sample of every condition in each batch (maybe with boolean flag for user to activate)
splimpDesignBatches <- function(sim, params, balancedBatches) {
  nConditions <- getParam(params, "nConditions")
  nBatches <- getParam(params, "nBatches")
  batch.size <- getParam(params, "batch.size")
  samples <- unique(colData(sim)$Sample)
  cell.names <- rownames(colData(sim))
  gene.names <- rownames(rowData(sim))
  if (nBatches == 1) {
    colData(sim)$Batch <- "Batch1"
  }
  else if (nBatches > 1){
    if (nConditions > 1 & balancedBatches & batch.size > 1) {
      batches <- rep(paste0("Batch", 1:nBatches), batch.size)
      conditions <- unique(colData(sim)$Condition)
      cond.assignment <- list()
      for (cond in conditions) {
        cond.assignment[[cond]] <- unique(colData(sim)[colData(sim)$Condition == cond,]$Sample)
      }
      batch.assignment <- list()
      # Add one sample from each condition to each batch in first iteration over batches
      for (batch in unique(batches)) {
        batch.assignment[[batch]] <- c()
        for (cond in conditions) {
          first.samp <- sample(cond.assignment[[cond]], 1)
          # Remove drawn sample from list
          cond.assignment[[cond]] <- cond.assignment[[cond]][!cond.assignment[[cond]] == first.samp]
          # Add drawn sample to batch
          batch.assignment[[batch]] <- c(batch.assignment[[batch]], first.samp)
        }
      }
      # Join list of both conditions
      cond.assignment <- Reduce(c, cond.assignment)
      cond.assignment.len <- length(cond.assignment)
      # Add remaining samples to the batches
      for (batch in unique(batches)) {
        draw <- sample(cond.assignment, floor(cond.assignment.len / length(unique(batches))), replace = FALSE)
        # Add to batches
        batch.assignment[[batch]] <- c(batch.assignment[[batch]], draw)
        # Removes drawn samples, s.t. samples can not be drawn twice
        cond.assignment <- cond.assignment[!cond.assignment %in% draw]
      }
      # Check whether all samples have been drawn
      # Not ready yet for unbalanced batches with differing amounts of samples
      stopifnot(length(cond.assignment) == 0)
      
      # Convert batch lists into data frame
      final.assignment <- data.frame(Sample = NULL, Batch = NULL)
      for (batch in unique(batches)) {
        final.assignment <- rbind(final.assignment, data.frame(Sample = batch.assignment[[batch]], Batch = batch))
      }
    } 
    else {
      stop("Currently only balanced batches available")
    } 

    colData(sim) <- merge(colData(sim), final.assignment, by = "Sample")
    # Cell names are lost during merging
    rownames(colData(sim)) <- cell.names
    
  }
  else {
    stop("At least one batch is necessary!")
  }
  return(sim)
}


# For each batch (whole simulation) and not for each sample separately (splatPopSimBatchEffects & splatSimBatchCellMeans)
# Pass global params with batch effects here (not params.intern)
splimpSimBatchEffects <- function(sim, params) {
  cell.names <- rownames(colData(sim))
  gene.names <- rownames(rowData(sim))
  nGenes <- length(gene.names)
  nBatches <- getParam(params, "nBatches")
  batch.facLoc <- getParam(params, "batch.facLoc")
  batch.facScale <- getParam(params, "batch.facScale")
  samples <- unique(colData(sim)$Sample)
  for (samp in samples) {
    # Get batch of sample
    batch <- unique(colData(sim)[, c('Sample', 'Batch')][colData(sim)$Sample == samp,])$Batch
    batch.num <- as.numeric(substring(batch, 6))
    # Generate batch effects
    batch.facs <- getLNormFactors(nGenes, 1, 0.5, batch.facLoc[batch.num], batch.facScale[batch.num])
    rowData(sim)[[paste(samp, "BatchFac", sep = "_")]] <- batch.facs
  }
  return(sim)
}

# Use Simulated_Means and BatchEffects to get batch.means.cell
splimpSimBatchCellMeans <- function(sim, params) {
  nBatches <- getParam(params, "nBatches")
  cell.names <- rownames(colData(sim))
  gene.names <- rownames(rowData(sim))
  samples <- unique(colData(sim)$Sample)
  batch.means.cell.list <- list()
  for (samp in samples) {
    gene.means <- metadata(sim)$Simulated_Means$Group1[, samp]
    nCells <- sum(colData(sim)$Sample == samp)
    if (nBatches > 1) {
      batch.facs.gene <- as.matrix(rowData(sim)[, paste(samp, "BatchFac", sep = "_")])
      batch.facs.cell <- as.matrix(batch.facs.gene[, rep(1, nCells)])
    } else {
      nGenes <- length(rownames(rowData(sim)))
      batch.facs.cell <- matrix(1, ncol = nCells, nrow = nGenes)
    }
    batch.means.cell.list[[samp]] <- batch.facs.cell * gene.means
  }
  batch.means.cell <- do.call(cbind, batch.means.cell.list)
  colnames(batch.means.cell) <- cell.names
  rownames(batch.means.cell) <- gene.names
  assays(sim)$BatchCellMeans <- batch.means.cell

  return(sim)
}

# Combine BatchCellMeans with LibSizeFactors to get BaseCellMeans (splatSimSingleCellMeans)
splimpSingleCellMeans <- function(sim) {
  nCells <- length(rownames(colData(sim)))
  cell.names <- rownames(colData(sim))
  gene.names <- rownames(rowData(sim))
  exp.lib.sizes <- colData(sim)$ExpLibSize
  batch.means.cell <- assays(sim)$BatchCellMeans
  cell.means.gene <- batch.means.cell

  # Do lib size normalization [results in colSums(cell.means.gene) == rep(1, X)]
  cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
  base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
  colnames(base.means.cell) <- cell.names
  rownames(base.means.cell) <- gene.names
  assays(sim)$BaseCellMeans <- base.means.cell

  return(sim)
}


# Combine BaseCellMeans with BCV to get CellMeans (splatSimBCVMeans)
splimpSimBCVMeans <- function(sim, params) {
  samples <- unique(colData(sim)$Sample)
  cell.names <- rownames(colData(sim))
  gene.names <- rownames(rowData(sim))
  nGenes <- length(gene.names)
  bcv.common <- getParam(params, "bcv.common")
  bcv.df <- getParam(params, "bcv.df")
  
  means.cell.list <- list()
  bcv.list <- list()
  for (samp in samples) {
    nCells <- sum(colData(sim)$Sample == samp)
    base.means.cell <- assays(sim[, sim$Sample == samp])$BaseCellMeans
    bcv.list[[samp]] <- (bcv.common + (1 / sqrt(base.means.cell))) * sqrt(bcv.df / rchisq(nGenes, df = bcv.df))
    means.cell.list[[samp]] <- matrix(rgamma(as.numeric(nGenes) * as.numeric(nCells),
                                              shape = 1 / (bcv.list[[samp]]^2),
                                              scale = base.means.cell * (bcv.list[[samp]]^2)),
                                       nrow = nGenes,
                                       ncol=nCells)
  }
  means.cell <- do.call(cbind, means.cell.list)
  bcv <- do.call(cbind, bcv.list)
  colnames(means.cell) <- cell.names
  rownames(means.cell) <- gene.names
  assays(sim)$BCV <- bcv
  assays(sim)$CellMeans <- means.cell

  return(sim)
}


# Simulate True counts from Poisson distribution (splatSimTrueCounts)
splimpSimTrueCounts <- function(sim) {
  cell.names <- rownames(colData(sim))
  gene.names <- rownames(rowData(sim))
  nCells <- length(cell.names)
  nGenes <- length(gene.names)
  cell.means <- assays(sim)$CellMeans

  true.counts <- matrix(rpois(as.numeric(nGenes) * as.numeric(nCells),
                              lambda = cell.means),
                        nrow = nGenes,
                        ncol = nCells)
  colnames(true.counts) <- cell.names
  rownames(true.counts) <- gene.names
  assays(sim)$TrueCounts <- true.counts
  assays(sim)$counts <- true.counts

  return(sim)
}


# Dropout is not used for this thesis
splimpSimDropout <- function(sim, params) {
  
  # Use function from splatter base
  sim <- splatSimDropout(sim, params)

  return(sim)
}


# Helper function to get LNormFactors
getLNormFactors <- function (n.facs, sel.prob, neg.prob, fac.loc, fac.scale)
{
  is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
  n.selected <- sum(is.selected)
  dir.selected <- (-1)^rbinom(n.selected, 1, neg.prob)
  facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
  # Reverse direction for factors that are less than one
  dir.selected[facs.selected < 1] <- -1 * dir.selected[facs.selected < 1]
  factors <- rep(1, n.facs)
  factors[is.selected] <- facs.selected^dir.selected
  return(factors)
}

# Get cli arguments
args <- commandArgs(trailingOnly = TRUE)
output <- args[1]
nGenes <- as.integer(args[2])
nCells <- as.integer(args[3])
seed <- as.integer(args[4])
preprocessing_threshold <- as.numeric(args[5])

set.seed(seed)



assignParams <- function() {
  vcf <- mockVCF(n.samples = 10)
  
  params <- newSplatPopParams(nGenes = nGenes,
                              eqtl.ES.shape = 0,
                              eqtl.ES.rate = 0,
                              out.prob = 0,
                              eqtl.n = 0,
                              eqtl.group.specific = 0,
                              eqtl.condition.specific = 0,
                              batch.size = 10,
                              batchCells = c(nCells / 10),
                              condition.prob = c(0.5, 0.5),
                              cde.prob = 0.025,
                              cde.downProb = 0.5,
                              cde.facLoc = 2.5,
                              cde.facScale = 0.4
  )
  return(list(params=params, vcf=vcf))
}

errorMAST <- function(sim) {
  # Do preprocessing according to the preprocessing process
  counts <- assay(sim)
  sim_filtered <- sim[rowSums(counts > 0) > (preprocessing_threshold * ncol(counts)), ]
  
  # Run MAST to check if it fails
  names(assays(sim))[1] <- "counts"
  
  exprsArray <- as.matrix(counts(sim))
  cData <- colData(sim)
  fData <- data.frame(Gene = rownames(exprsArray))
  
  exprsArray.log <- log2(exprsArray + 1)
  
  sca <- FromMatrix(exprsArray = exprsArray.log, cData = cData, fData = fData, check_sanity = F)
  
  colData(sca)$Batch <- factor(colData(sca)$Batch)
  colData(sca)$Condition <- factor(colData(sca)$Condition)
  colData(sca)$Sample <- factor(colData(sca)$Sample)
  
  tryCatch({
    zlmCond <- zlm(~ Condition + (1 | Sample), sca, method = "glmer", ebayes = F, strictConvergence = F)
    print("No errors with MAST")
    return (FALSE)
  }, error = function(e) {
    return (TRUE)
  })
}


assignedParams <- assignParams()
sim <- splimpSimulate(vcf = assignedParams$vcf, params = assignedParams$params)

while (errorMAST(sim)) {
  assignedParams <- assignParams()
  sim <- splimpSimulate(vcf = assignedParams$vcf, params = assignedParams$params)
}

assays(sim)$BCV <- NULL
assays(sim)$BaseCellMeans <- NULL
assays(sim)$BatchCellMeans <- NULL
assays(sim)$CellMeans <- NULL
assays(sim)$TrueCounts <- NULL

write_h5ad(sim, output)
