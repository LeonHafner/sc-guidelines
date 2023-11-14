#!/bin/bash

#SBATCH --job-name=pvalues_dataset
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=1-20


python3 tools/pvalues.py \
	--input data \
	--output data/12_pvalues/pvalues_${SLURM_ARRAY_TASK_ID}.tsv \
	--run ${SLURM_ARRAY_TASK_ID}
