#!/bin/bash

#SBATCH --job-name=distinct_kang2018
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err

Rscript tools/distinct.R \
       --input real-data/data/2_preprocessed/preprocessed.h5ad \
       --output real-data/data/6_distinct/distinct.tsv \
       --scenario kang2018
