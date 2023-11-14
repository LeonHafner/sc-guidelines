#!/bin/bash

#SBATCH --job-name=atlas_scvi
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


python3 tools/scvi-de.py \
	--input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/10_scvi/scvi_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario atlas
