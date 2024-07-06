#!/usr/bin/env python3

import scanpy as sc
import anndata as ad
import random

adata = sc.read("${input_anndata}")

adatas = []
for batch in adata.obs['Batch'].unique():
    adata_batch = adata[adata.obs['Batch'] == batch]
    condition = adata_batch.obs['Condition'].tolist()
    random.shuffle(condition)
    adata_batch.obs['Condition'] = condition
    adatas.append(adata_batch)

adata = ad.concat(adatas)
adata.obs['Sample'] = adata.obs['Batch'].astype(str) + '-' + adata.obs['Condition'].astype(str)

adata.write('kang2018.h5ad')