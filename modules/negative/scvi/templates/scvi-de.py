#!/usr/bin/env python3

import pandas as pd
from math import ceil
import numpy as np
import scanpy as sc
import scvi
from tqdm import tqdm

scvi.settings.seed = 0

file_in = "${input_anndata}"
scenario = "${meta.scenario}"


print(f'Input: {file_in}')
print(f'Scenario: {scenario}')

adata = sc.read(file_in)

if scenario in ["atlas", "atlas_hvg", "atlas-less-de"]:
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Batch", "Sample"],
    )
elif scenario in ["dataset", "dataset_hvg", "dataset-less-de"]:
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Sample"],
    )
elif scenario in ["atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de"]:
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Batch", "Sample"],
    )
elif scenario == "atlas-negative":
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Batch", "Sample"],
    )
elif scenario == "kang2018":
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Batch", "Sample"],
    )
elif scenario in ["dataset-ub-cells", "dataset-ub-cells_hvg", "dataset-ub-cells-less-de"]:
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Sample"],
    )

model = scvi.model.SCVI(adata)


# Update batch size based on this issue: https://github.com/scverse/scvi-tools/issues/2214
batch_size = 128
while ceil(adata.n_obs * 0.9) % batch_size == 1:
    batch_size += 1
    print(f'Updating batch size to {batch_size}')
print(f'Using batch size of {batch_size}')

model.train(batch_size=batch_size)

group1 = "Condition1"
group2 = "Condition2"

for cutoff in tqdm(np.linspace(0, 1, 2001)):
    de = model.differential_expression(
        groupby="Condition", group1=group1, group2=group2, fdr_target=round(cutoff, 4)
    )
    de.sort_index().to_csv(f'result-{str(round(cutoff, 4)).replace(".", "_")}.tsv', sep = "\\t")

