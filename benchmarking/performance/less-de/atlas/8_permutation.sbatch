#!/bin/bash

#SBATCH --job-name=atlas_permutation
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


python3 tools/permutation.py \
	--input data/3_pseudobulked/pseudobulked_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/8_permutation/permutation_${SLURM_ARRAY_TASK_ID}.tsv \
	-n 10000 \
	-n_max 100000
