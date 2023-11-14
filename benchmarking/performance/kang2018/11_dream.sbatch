#!/bin/bash

#SBATCH --job-name=dream_kang2018
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


Rscript tools/dream.R \
        --input data/3_pseudobulked/pseudobulked.h5ad \
        --output data/11_dream/dream.tsv \
        --scenario kang2018 \
        --threads 8
