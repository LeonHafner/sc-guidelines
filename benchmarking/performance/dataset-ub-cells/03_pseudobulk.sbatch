#!/bin/bash

#SBATCH --job-name=pseudobulk_dataset-ub-cells
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


python3 tools/pseudobulk.py \
       --input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
       --output data/3_pseudobulked/pseudobulked_${SLURM_ARRAY_TASK_ID}.h5ad \
       --groupby Sample Condition
