#!/bin/bash

#SBATCH --job-name=runtime_benchmark
#SBATCH --cpus-per-task=8 
#SBATCH --mem=8gb 
#SBATCH -t 4-00:00:00 
#SBATCH -o out_slurm/%A_%a.out
#SBATCH -e out_slurm/%A_%a.err
#SBATCH --array=0-27


python3 /nfs/home/students/l.hafner/fopra/runtime/runtime_benchmark.py ${1} ${SLURM_ARRAY_TASK_ID}
