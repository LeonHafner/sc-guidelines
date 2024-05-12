#!/usr/bin/env python3
   
import scanpy as sc

ad = sc.read_h5ad("${input_anndata}")

sc.pp.highly_variable_genes(
    adata=ad,
    n_top_genes=int(ad.shape[1] * float("${hvg_ratio}")),
    flavor='seurat_v3',
    subset=True,
)

ad.write_h5ad("${meta.scenario}_${meta.run}.h5ad")