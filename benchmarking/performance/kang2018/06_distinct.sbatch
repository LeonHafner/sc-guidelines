#!/bin/bash

#SBATCH --job-name=kang2018_distinct
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


Rscript tools/distinct.R \
       --input data/2_preprocessed/preprocessed.h5ad \
       --output data/6_distinct/distinct.tsv \
       --scenario kang2018
