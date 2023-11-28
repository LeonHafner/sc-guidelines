#!/bin/bash

#SBATCH --job-name=kang2018_pseudobulk
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb 
#SBATCH -t 00:02:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/pseudobulk.py \
       --input data/2_preprocessed/preprocessed.h5ad \
       --output data/3_pseudobulked/pseudobulked.h5ad \
       --groupby ind sample stim
