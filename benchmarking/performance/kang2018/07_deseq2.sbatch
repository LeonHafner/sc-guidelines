#!/bin/bash

#SBATCH --job-name=kang2018_deseq2
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err

Rscript tools/deseq2.R \
        --input data/3_pseudobulked/pseudobulked.h5ad \
        --output data/7_deseq2/deseq2.tsv \
        --scenario kang2018
