#!/usr/bin/env python
import argparse
import pandas as pd
import scanpy as sc
import scvi


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

if scenario == "kang2018":
    adata.obs.columns = ["Batch", "Condition", "cluster", "cell", "multiplets", "Sample"]

if scenario == "atlas":
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Batch", "Sample"],
    )
elif scenario == "dataset":
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Sample"],
    )
elif scenario == "atlas-ub-conditions":
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
elif scenario == "dataset-ub-cells":
    scvi.model.SCVI.setup_anndata(
        adata,
        categorical_covariate_keys=["Sample"],
    )

model = scvi.model.SCVI(adata)
model.train()

group1 = "Condition1"
group2 = "Condition2"

if scenario == "kang2018":
    group1 = "stim"
    group2 = "ctrl"

de = model.differential_expression(
    groupby="Condition", group1=group1, group2=group2
)

de.sort_index().to_csv(file_out, sep = "\t")

