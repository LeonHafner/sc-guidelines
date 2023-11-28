#!/bin/bash

#SBATCH --job-name=kang2018_preprocessing
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb 
#SBATCH -t 02:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/preprocessing.py \
        --input data/1_raw/raw.h5ad \
        --output data/2_preprocessed/preprocessed.h5ad \
        --threshold 0.1
