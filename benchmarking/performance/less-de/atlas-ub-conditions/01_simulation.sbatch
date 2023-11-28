#!/bin/bash

#SBATCH --job-name=sim_atlas-ub-conditions
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=20


Rscript tools/simulation.R \
	--output data/1_raw/raw_${SLURM_ARRAY_TASK_ID}.h5ad \
	--scenario atlas-ub-conditions_less_de
