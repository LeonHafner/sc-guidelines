#!/bin/bash

#SBATCH --job-name=pb_fixed_pvalues_dataset-ub-cells
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


python3 tools/pvalues.py \
	--input data \
	--output data/12_pvalues/pvalues_${SLURM_ARRAY_TASK_ID}.tsv \
	--run ${SLURM_ARRAY_TASK_ID}
