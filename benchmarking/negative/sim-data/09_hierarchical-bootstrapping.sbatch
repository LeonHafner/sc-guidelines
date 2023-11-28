#!/bin/bash

#SBATCH --job-name=hierarchical_bootstrapping_atlas-negative
#SBATCH -p sesame_street
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/hierarchical_bootstrapping.py \
	--input data/2_preprocessed/preprocessed.h5ad \
	--output1 data/9_hierarchical-bootstrapping/hierarchical-bootstrapping.tsv \
	--output2 data/9_hierarchical-bootstrapping/hierarchical-bootstrapping_adv.tsv \
	--scenario atlas-negative
