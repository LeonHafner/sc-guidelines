#!/usr/bin/env python
# coding: utf-8
import sys
import time
import itertools
from tqdm import tqdm
import subprocess
import pandas as pd

run = int(sys.argv[1])
mode = 'genes_fixed'

if mode == 'cells_fixed':
    n_genes = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
               1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000,
               5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000]
    n_cells = [1000]
elif mode == 'genes_fixed':
    n_cells = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
               1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000,
               5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000]
    n_genes = [1000]
else:
    raise Error("Wrong argument")

combination = list(itertools.product(n_genes, n_cells))[int(sys.argv[2])]

data = {
    'index': [],
    'sim_start': [],
    'sim_end': [],
    'preprocess_start': [],
    'preprocess_end': [],
    'pseudobulk_start': [],
    'pseudobulk_end': [],
    'de-genes_start': [],
    'de-genes_end': [],
    'mast_start': [],
    'mast_end': [],
    'distinct_start': [],
    'distinct_end': [],
    'deseq2_start': [],
    'deseq2_end': [],
    'permutation_start': [],
    'permutation_end': [],
    'hierarchical-bootstrapping_start': [],
    'hierarchical-bootstrapping_end': [],
    'scvi_start': [],
    'scvi_end': [],
    'dream_start': [],
    'dream_end': [],
}

n_g, n_c = combination
filename = f'sim_{n_g}_{n_c}_run{run}'
basepath = 'path/to/base-dir'
data['index'].append(filename)

# Simulation
data['sim_start'].append(time.time())
subprocess.call(["conda", "run","-n", "<environment>", "-v",
                 "Rscript", "simulation-runtime.R",
                 f'{filename}.h5ad', str(n_g), str(n_c), mode])
data['sim_end'].append(time.time())


# Preprocessing (reducing genes)
print("Preprocessing\n")
data['preprocess_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'python3', f'{basepath}/tools/preprocessing.py',
                 '--input', f'{basepath}/runtime/data/{mode}/1_sims/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/2_preprocessed/{filename}.h5ad'])
data['preprocess_end'].append(time.time())


# Pseudobulking
print("Pseudobulking\n")
data['pseudobulk_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'python3', f'{basepath}/tools/pseudobulk.py',
                 '--input', f'{basepath}/runtime/data/{mode}/2_preprocessed/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/3_pseudobulked/{filename}.h5ad'])
data['pseudobulk_end'].append(time.time())


# Get DE genes
print("DE genes\n")
data['de-genes_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'python3', f'{basepath}/tools/de-genes.py',
                 '--input', f'{basepath}/runtime/data/{mode}/2_preprocessed/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/4_de-genes/{filename}.tsv'])
data['de-genes_end'].append(time.time())


# MAST
print("MAST\n")
data['mast_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'Rscript', f'{basepath}/tools/mast.R',
                 '--input', f'{basepath}/runtime/data/{mode}/2_preprocessed/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/5_mast/{filename}.tsv',
                 '--scenario', 'dataset'])
data['mast_end'].append(time.time())


# distinct
print("distinct\n")
data['distinct_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'Rscript', f'{basepath}/tools/distinct.R',
                 '--input', f'{basepath}/runtime/data/{mode}/2_preprocessed/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/6_distinct/{filename}.tsv',
                 '--scenario', 'dataset'])
data['distinct_end'].append(time.time())


# DESeq2
print("DESeq2\n")
data['deseq2_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'Rscript', f'{basepath}/tools/deseq2.R',
                 '--input', f'{basepath}/runtime/data/{mode}/3_pseudobulked/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/7_deseq2/{filename}.tsv',
                 '--scenario', 'dataset'])
data['deseq2_end'].append(time.time())


# Permutation
print("Permutation\n")
data['permutation_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'python3', f'{basepath}/tools/permutation.py',
                 '--input', f'{basepath}/runtime/data/{mode}/3_pseudobulked/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/8_permutation/{filename}.tsv'])
data['permutation_end'].append(time.time())


# Hierarchical bootstrapping
print("Hierarchical Bootstrapping\n")
data['hierarchical-bootstrapping_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'python3', f'{basepath}/tools/hierarchical_bootstrapping.py',
                 '--input', f'{basepath}/runtime/data/{mode}/2_preprocessed/{filename}.h5ad',
                 '--output1', f'{basepath}/runtime/data/{mode}/9_hierarchical-bootstrapping/{filename}.tsv',
                 '--output2', f'{basepath}/runtime/data/{mode}/9_hierarchical-bootstrapping/{filename}_adv.tsv',
                 '--scenario', 'dataset'])
data['hierarchical-bootstrapping_end'].append(time.time())
    


# scVI
print("scVI\n")
data['scvi_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', 'scvi_test', '-v',
                 'python3', f'{basepath}/tools/scvi-de.py',
                 '--input', f'{basepath}/runtime/data/{mode}/2_preprocessed/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/10_scvi/{filename}.tsv',
                 '--scenario', 'dataset'])
data['scvi_end'].append(time.time())


# DREAM
print("DREAM\n")
data['dream_start'].append(time.time())
subprocess.call(['conda', 'run', '-n', '<environment>', '-v',
                 'Rscript', f'{basepath}/tools/dream.R',
                 '--input', f'{basepath}/runtime/data/{mode}/3_pseudobulked/{filename}.h5ad',
                 '--output', f'{basepath}/runtime/data/{mode}/11_dream/{filename}.tsv',
                 '--scenario', 'dataset',
                 '--threads', '8'])
data['dream_end'].append(time.time())


data = pd.DataFrame.from_dict(data)


data.to_csv(f'{basepath}/runtime/data/{mode}/results/{filename}.tsv', sep = '\t', index = False)

