#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import itertools
import random
from tqdm import tqdm
from scipy.stats import spearmanr

random.seed(int("${meta.run}"))

path_anndata = "${anndata}"
run = "${meta.run}"

adata = sc.read(path_anndata)

result_intra = []
for sample in tqdm(adata.obs['Sample'].unique()):
    adata_sample = adata[adata.obs['Sample'] == sample]
    print(sample, adata_sample.shape[0])
    for i, j in itertools.combinations(range(adata_sample.n_obs), r=2):
        result_intra.append(('intra', run, spearmanr(adata_sample.X[i], adata_sample.X[j]).statistic))

result_inter = []
for _ in tqdm(range(1000)):
    cells = []
    for sample in adata.obs['Sample'].unique():
        adata_sample = adata[adata.obs['Sample'] == sample]
        cells.append(adata_sample.X[random.randint(0, adata_sample.n_obs - 1), :])
    for cell_1, cell_2 in itertools.combinations(cells, r=2):
        result_inter.append(('inter', run, spearmanr(cell_1, cell_2).statistic))


result = result_intra + result_inter
pd.DataFrame(result, columns=['type', 'run', 'correlation']).to_csv("correlations.tsv", sep="\\t", index=False)