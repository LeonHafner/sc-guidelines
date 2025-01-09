#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd

parser = argparse.ArgumentParser(description="Perform DEA with scanpy's rank_genes_group()")
parser.add_argument("--input", required=True, help="Path to the input AnnData file.")
parser.add_argument("--output", required=True, help="Path to the output CSV file.")
args = parser.parse_args()

adata = sc.read(args.input)

sc.pp.log1p(adata)
sc.tl.rank_genes_groups(adata, groupby="Condition", method="t-test")

results = adata.uns["rank_genes_groups"]
group = "Condition1" if "Condition1" in [name for name, _ in results["names"].dtype.descr] else "lung_adenocarcinoma"

genes = results['names'][group]
pvals = results['pvals'][group]

gene_to_pval = {gene: pval for gene, pval in zip(genes, pvals)}

pd.DataFrame({
    "gene": adata.var_names,
    "pvalue": [gene_to_pval[gene] if gene in gene_to_pval else pd.NA for gene in adata.var_names],
}).to_csv(args.output, index=False, sep='\t')
