#!/bin/bash

#SBATCH --job-name=kang2018_permutation
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/permutation.py \
        --input data/3_pseudobulked/pseudobulked.h5ad \
        --output data/8_permutation/permutation.tsv \
        -n 10000 \
        -n_max 100000
