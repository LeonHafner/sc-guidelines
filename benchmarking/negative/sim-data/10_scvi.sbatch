#!/bin/bash

#SBATCH --job-name=scvi_atlas-negative
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/scvi-de.py \
	--input data/2_preprocessed/preprocessed.h5ad \
	--output data/10_scvi/scvi.tsv \
	--scenario atlas-negative
