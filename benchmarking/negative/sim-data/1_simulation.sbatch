#!/bin/bash

#SBATCH --job-name=sim_atlas_negative
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


Rscript tools/simulation.R \
	--output data/1_raw/raw.h5ad \
	--scenario atlas-negative
