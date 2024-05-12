#!/usr/bin/env python3

import scanpy as sc
import numpy as np


file_in = "${input_anndata}"
file_out = "de_${meta.scenario}_${meta.run}.tsv"
scenario = "${meta.scenario}"
reactome = "${reactome}"

print(f'Input: {file_in}')
print(f'Output: {file_out}')
print(f'Scenario: {scenario}')
print(f'Reactome: {reactome}')

adata = sc.read_h5ad(file_in)

if scenario == 'kang2018':
    print(f'Reactome: {reactome}')
    # Load reactome genes
    with open(reactome, "r") as f:
        de_genes = f.read().split("\\n")
    while "" in de_genes:
        de_genes.remove("")
    
    # Update de-genes in adata.var
    adata.var['ConditionDE.Condition1'] = 1.0
    adata.var['ConditionDE.Condition2'] = 1.0
    mask = np.isin(adata.var.index.tolist(), de_genes)

    # Set de-genes to arbitrary value other than 1.0
    adata.var.loc[mask, 'ConditionDE.Condition1'] = 10.0

adata.var[['ConditionDE.Condition1', 'ConditionDE.Condition2']].to_csv(file_out, sep = '\\t')

