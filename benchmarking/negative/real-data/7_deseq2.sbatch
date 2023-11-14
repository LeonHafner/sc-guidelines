#!/bin/bash

#SBATCH --job-name=deseq2_kang2018
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err

Rscript tools/deseq2.R \
	--input real-data/data/3_pseudobulked/pseudobulked.h5ad \
	--output real-data/data/7_deseq2/deseq2.tsv \
	--scenario kang2018
