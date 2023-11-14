#!/bin/bash

#SBATCH --job-name=sim_dataset-ub-cells
#SBATCH --cpus-per-task=16
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=1-20


Rscript tools/simulation.R \
	--output data/1_raw/raw_${SLURM_ARRAY_TASK_ID}.h5ad \
	--scenario dataset-ub-cells_less_de
