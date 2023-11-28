#!/bin/bash

#SBATCH --job-name=preprocessing_kang2018
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/preprocessing.py \
	--input real-data/data/1_raw/raw.h5ad \
	--output real-data/data/2_preprocessed/preprocessed.h5ad \
	--threshold 0.1
