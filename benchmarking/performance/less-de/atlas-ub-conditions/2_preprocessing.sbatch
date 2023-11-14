#!/bin/bash

#SBATCH --job-name=preprocessing_atlas-ub-conditions
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=20


python3 tools/unbalance-atlas.py \
	--input data/1_raw/raw_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad


python3 tools/preprocessing.py \
	--input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--threshold 0.1
