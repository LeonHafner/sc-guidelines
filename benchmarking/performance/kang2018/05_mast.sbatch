#!/bin/bash

#SBATCH --job-name=kang2018_mast
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 04:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


Rscript tools/mast.R \
        --input data/2_preprocessed/preprocessed.h5ad \
        --output data/5_mast/mast.tsv \
        --scenario kang2018
