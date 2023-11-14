#!/bin/bash

#SBATCH --job-name=precision-recall_dataset-ub-cells
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 14-0:00 
#SBATCH -o out_slurm/%j.out
#SBATCH -e out_slurm/%j.err
#SBATCH --array=1-20


python3 tools/precision-recall.py \
	--pvalues data/12_pvalues/pvalues_${SLURM_ARRAY_TASK_ID}.tsv \
	--de_genes data/4_de-genes/de-genes_${SLURM_ARRAY_TASK_ID}.tsv \
	--output data/13_precision-recall \
	--run ${SLURM_ARRAY_TASK_ID}
