#!/bin/bash

#SBATCH --job-name=distinct_dataset-ub-cells
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


Rscript tools/distinct.R \
       --input performance_less_de/dataset-ub-cells/data/2_preprocessed/preprocessed_${SLURM_ARRAY_TASK_ID}.h5ad \
       --output performance_less_de/dataset-ub-cells/data/6_distinct/distinct_${SLURM_ARRAY_TASK_ID}.tsv \
       --scenario dataset-ub-cells
