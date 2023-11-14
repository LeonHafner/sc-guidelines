#!/bin/bash

#SBATCH --job-name=mast_dataset
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


Rscript tools/mast.R \
	--input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/5_mast/mast_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario dataset
