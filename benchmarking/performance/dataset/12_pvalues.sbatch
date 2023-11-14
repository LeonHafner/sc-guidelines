#!/bin/bash

#SBATCH --job-name=pvalues_dataset
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


python3 tools/pvalues.py \
	--input data \
	--output data/12_pvalues/pvalues_${SLURM_ARRAY_TASK_ID}.tsv \
	--run ${SLURM_ARRAY_TASK_ID}
