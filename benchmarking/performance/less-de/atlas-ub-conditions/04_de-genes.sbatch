#!/bin/bash

#SBATCH --job-name=de-genes_atlas-ub-conditions
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=12,20


python3 tools/de-genes.py \
       --input data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
       --output data/4_de-genes/de-genes_${SLURM_ARRAY_TASK_ID}.tsv
