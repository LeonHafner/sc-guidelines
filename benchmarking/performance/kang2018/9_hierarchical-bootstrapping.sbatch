#!/bin/bash

#SBATCH --job-name=hierarchical_bootstrapping_kang2018
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err


python3 tools/hierarchical_bootstrapping.py \
        --input data/2_preprocessed/preprocessed.h5ad \
        --output1 data/9_hierarchical-bootstrapping/hierarchical-bootstrapping.tsv \
        --output2 data/9_hierarchical-bootstrapping/hierarchical-bootstrapping_adv.tsv \
        --scenario kang2018
