#!/bin/bash

#SBATCH --job-name=mast_dataset-ub-cells
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=1-20


Rscript tools/mast.R \
	--input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/5_mast/mast_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario dataset-ub-cells
