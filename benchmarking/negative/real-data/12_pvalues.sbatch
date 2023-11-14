#!/bin/bash

#SBATCH --job-name=pvalues_kang2018
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/pvalues.py \
	--input real-data/data \
	--output real-data/data/12_pvalues/pvalues.tsv \
	--run 1
