#!/usr/bin/env python3

import scanpy as sc
from scipy.stats import spearmanr
from tqdm import tqdm
import random

adata = sc.read("${anndata}")
adata_final = adata.copy()

correlations = abs(spearmanr(adata.X[:, 0], adata.X[:, 1:]).statistic)

uncorrelated_genes = []
for _ in tqdm(range(500)):
    gene_name = random.sample(adata.var_names.tolist(), k=1)[0]
    uncorrelated_genes.append(gene_name)
    gene_index = adata.var_names.tolist().index(gene_name)
    
    correlation = correlations[gene_index, :]

    adata = adata[:, correlation < 0.1]
    correlations = correlations[correlation < 0.1, :][:, correlation < 0.1]
    if adata.shape[1] == 0:
        break

adata = adata_final[:, uncorrelated_genes]
adata.write_h5ad("uncorrelated_genes.h5ad")