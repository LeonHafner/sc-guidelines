#!/bin/bash

#SBATCH --job-name=dream_atlas-ub-conditions
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=12,20


Rscript tools/dream.R \
	--input data/3_pseudobulked/pseudobulked_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/11_dream/dream_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario atlas-ub-conditions \
	--threads 8
