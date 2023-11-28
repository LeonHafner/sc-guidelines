#!/bin/bash

#SBATCH --job-name=permutation_atlas-negative
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/permutation.py \
	--input data/3_pseudobulked/pseudobulked.h5ad \
	--output data/8_permutation/permutation.tsv \
	-n 10000 \
	-n_max 100000
