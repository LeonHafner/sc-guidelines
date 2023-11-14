#!/bin/bash

#SBATCH --job-name=hierarchical_bootstrapping_dataset
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


python3 tools/hierarchical_bootstrapping.py \
	--input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output1 data/9_hierarchical-bootstrapping/hierarchical-bootstrapping_${SLURM_ARRAY_TASK_ID}.tsv \
	--output2 data/9_hierarchical-bootstrapping/hierarchical-bootstrapping_adv_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario dataset
