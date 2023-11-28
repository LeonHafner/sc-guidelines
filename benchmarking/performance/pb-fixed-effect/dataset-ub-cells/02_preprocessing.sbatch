#!/bin/bash

#SBATCH --job-name=pb_fixed_preprocessing_dataset-ub-cells
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


python3 tools/preprocessing.py \
	--input data/1_raw/raw_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--threshold 0.1
