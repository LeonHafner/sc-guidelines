#!/bin/bash

#SBATCH --job-name=mast_atlas-ub-conditions
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=12,20


Rscript tools/mast.R \
	--input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/5_mast/mast_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario atlas-ub-conditions
