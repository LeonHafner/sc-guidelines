#!/usr/bin/env python3

import argparse
from math import ceil
import pandas as pd
import scanpy as sc
import scvi

scvi.settings.seed = 0

parser = argparse.ArgumentParser(
                    prog='scVI approach',
                    description='scVI differential expression testing approach')

parser.add_argument('--input', help='input h5ad', required=True)
parser.add_argument('--output', help='output tsv', required=True)
parser.add_argument('--scenario', help = "Scenario of data", required=True)


args = parser.parse_args()
file_in = args.input
file_out = args.output
scenario = args.scenario


print(f'Input: {file_in}')
print(f'Output: {file_out}')
print(f'Scenario: {scenario}')

adata = sc.read(file_in)

if scenario in ["atlas", "atlas_hvg", "atlas-less-de"]:
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key="Batch",
        categorical_covariate_keys=["Sample"],
    )
elif scenario in ["dataset", "dataset_hvg", "dataset-less-de"]:
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Sample"],
    )
elif scenario in ["atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de"]:
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key="Batch",
        categorical_covariate_keys=["Sample"],
    )
elif scenario == "atlas-negative":
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key="Batch",
        categorical_covariate_keys=["Sample"],
    )
elif scenario == "kang2018":
    scvi.model.SCVI.setup_anndata(
        adata,
        batch_key="Batch",
        categorical_covariate_keys=["Sample"],
    )
elif scenario in ["dataset-ub-cells", "dataset-ub-cells_hvg", "dataset-ub-cells-less-de"]:
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Sample"],
    )
elif scenario == "luca":
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

if scenario == "luca":
    group1 = "lung_adenocarcinoma"
    group2 = "squamous_cell_lung_carcinoma"
else:
    group1 = "Condition1"
    group2 = "Condition2"

de = model.differential_expression(
    groupby="Condition",
    group1=group1,
    group2=group2,
    batch_correction=True,
)

de.sort_index().to_csv(file_out, sep = "\t")

