#!/bin/bash

#SBATCH --job-name=deseq2_atlas-negative
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


Rscript tools/deseq2.R \
	--input data/3_pseudobulked/pseudobulked.h5ad \
	--output data/7_deseq2/deseq2.tsv \
	--scenario atlas-negative
