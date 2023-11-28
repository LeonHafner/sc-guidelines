#!/bin/bash

#SBATCH --job-name=deseq2_dataset
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=1-20


Rscript tools/deseq2.R \
	--input data/3_pseudobulked/pseudobulked_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/7_deseq2/deseq2_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario dataset
