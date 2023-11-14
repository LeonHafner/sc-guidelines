#!/bin/bash

#SBATCH --job-name=scvi_kang2018
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/scvi-de.py \
        --input data/2_preprocessed/preprocessed.h5ad \
        --output data/10_scvi/scvi.tsv \
        --scenario kang2018
