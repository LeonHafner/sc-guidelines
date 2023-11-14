#!/bin/bash

#SBATCH --job-name=hierarchical_bootstrapping_atlas-ub-conditions
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=12,20


python3 tools/hierarchical_bootstrapping.py \
	--input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
	--output1 data/9_hierarchical-bootstrapping/hierarchical-bootstrapping_${SLURM_ARRAY_TASK_ID}.tsv \
	--output2 data/9_hierarchical-bootstrapping/hierarchical-bootstrapping_adv_${SLURM_ARRAY_TASK_ID}.tsv \
	--scenario atlas-ub-conditions
