#!/usr/bin/env python3

import argparse
import scanpy as sc
from typing import Callable, Sequence, Union
from anndata import AnnData
import numpy as np
import pandas as pd
from operator import and_
from functools import reduce


parser = argparse.ArgumentParser(
                    prog='Pseudobulking Step',
                    description='Program applies pseudobulk function to the input file and outputs a h5ad ')

parser.add_argument('--input', help='input h5ad', required=True)
parser.add_argument('--output', help='output h5ad', required=True)
parser.add_argument('--scenario', '-s', help='scenario of the data', required=True)

args = parser.parse_args()
file_in = args.input
file_out = args.output
scenario = args.scenario

print(f'Input: {file_in}')
print(f'Output: {file_out}')
print(f'Scenario: {scenario}')


def pseudobulk(
    adata: AnnData,
    *,
    groupby: Union[str, Sequence[str]],
    aggr_fun: Callable = np.sum,
    min_obs: int = 10,
) -> AnnData:
    """
    Calculate Pseudobulk of groups

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        One or multiple columns to group by
    aggr_fun
        Callback function to calculate pseudobulk. Must be a numpy ufunc supporting
        the `axis` attribute.
    min_obs
        Exclude groups with less than `min_obs` observations

    Returns
    -------
    New anndata object with same vars as input, but reduced number of obs.
    """
    if isinstance(groupby, str):
        groupby = [groupby]

    combinations = adata.obs.loc[:, groupby].drop_duplicates()

    # precompute masks
    masks = {}
    for col in groupby:
        masks[col] = {}
        for val in combinations[col].unique():
            masks[col][val] = adata.obs[col] == val

    expr_agg = []
    obs = []

    for comb in combinations.itertuples(index=False):
        mask = reduce(and_, (masks[col][val] for col, val in zip(groupby, comb)))
        if np.sum(mask) < min_obs:
            continue
        expr_row = aggr_fun(adata.X[mask, :], axis=0)
        obs_row = comb._asdict()
        obs_row["n_obs"] = np.sum(mask)
        # convert matrix to array if required (happens when aggregating spares matrix)
        try:
            expr_row = expr_row.A1
        except AttributeError:
            pass
        obs.append(obs_row)
        expr_agg.append(expr_row)

    return AnnData(
        X=np.vstack(expr_agg),
        var=adata.var,
        obs=pd.DataFrame.from_records(obs),
    )



adata = sc.read_h5ad(file_in)

match scenario:
    case 'atlas' | 'atlas_hvg' | 'atlas-less-de' | 'atlas-negative':
        groupby = ['Batch', 'Sample', 'Condition']
    case 'dataset' | 'dataset_hvg' | 'dataset-less-de' | 'luca':
        groupby = ['Sample', 'Condition']
    case 'atlas-ub-conditions' | 'atlas-ub-conditions_hvg' | 'atlas-ub-conditions-less-de':
        groupby = ['Batch', 'Sample', 'Condition']
    case 'dataset-ub-cells' | 'dataset-ub-cells_hvg' | 'dataset-ub-cells-less-de':
        groupby = ['Sample', 'Condition']
    case 'kang2018':
        groupby = ['Batch', 'Sample', 'Condition']
    case _:
        raise ValueError("Scenario not defined!")

adata_pb = pseudobulk(adata = adata,
                      groupby = groupby,
                      aggr_fun = np.sum,
                      min_obs = 0)

adata_pb.write_h5ad(file_out)
