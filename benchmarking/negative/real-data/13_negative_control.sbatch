#!/bin/bash

#SBATCH --job-name=negative-control_atlas-negative
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/negative-control.py \
	--pvalues real-data/data/12_pvalues/pvalues.tsv \
	--output real-data/data/13_negative-control
