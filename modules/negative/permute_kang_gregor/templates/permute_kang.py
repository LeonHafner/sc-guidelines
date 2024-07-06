#!/usr/bin/env python3

from math import ceil, floor
import scanpy as sc
import anndata as ad
import numpy as np
import random

adata = sc.read("${input_anndata}")

adatas = []
for batch in adata.obs['Batch'].unique():
    adata_batch = adata[adata.obs['Batch'] == batch].copy()
    cond1_idx = np.where(adata_batch.obs['Condition'] == "Condition1")[0].tolist()
    cond2_idx = np.where(adata_batch.obs['Condition'] == "Condition2")[0].tolist()
    random.shuffle(cond1_idx)
    random.shuffle(cond2_idx)
    condition_new = np.empty((adata_batch.shape[0],), dtype=object)
    condition_new[cond1_idx[:ceil(len(cond1_idx) / 2)] + cond2_idx[:ceil(len(cond2_idx) / 2)]] = "Condition1"
    condition_new[cond1_idx[floor(len(cond1_idx) / 2):] + cond2_idx[floor(len(cond2_idx) / 2):]] = "Condition2"
    adata_batch.obs['Condition'] = condition_new
    adatas.append(adata_batch)

adata = ad.concat(adatas)
adata.obs['Sample'] = adata.obs['Batch'].astype(str) + '-' + adata.obs['Condition'].astype(str)

adata.write('kang2018.h5ad')