#!/bin/bash

#SBATCH --job-name=atlas_dream
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


Rscript tools/dream.R \
	--input data/3_pseudobulked/pseudobulked_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/11_dream/dream_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario atlas \
	--threads 8
