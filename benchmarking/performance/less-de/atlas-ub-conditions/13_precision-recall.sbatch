#!/bin/bash

#SBATCH --job-name=precision-recall_atlas-ub-conditions
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=1-20


python3 tools/precision-recall.py \
	--pvalues data/12_pvalues/pvalues_${SLURM_ARRAY_TASK_ID}.tsv \
	--de_genes data/4_de-genes/de-genes_${SLURM_ARRAY_TASK_ID}.tsv \
	--output data/13_precision-recall \
	--run ${SLURM_ARRAY_TASK_ID}
