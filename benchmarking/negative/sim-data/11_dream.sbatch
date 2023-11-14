#!/bin/bash

#SBATCH --job-name=dream_atlas-negative
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


Rscript tools/dream.R \
	--input data/3_pseudobulked/pseudobulked.h5ad \
	--output data/11_dream/dream.tsv \
	--scenario atlas-negative \
	--threads 16
