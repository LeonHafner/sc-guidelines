#!/usr/bin/env python3

import argparse
from anndata import AnnData
import scanpy as sc
import numpy as np
from typing import Sequence, Tuple, Callable
import pandas as pd
from tqdm import tqdm

np.random.seed(0)

parser = argparse.ArgumentParser(
                    prog='Hierarchical Bootstrapping',
                    description='Conducts Hierarchical Bootstrapping on an anndata object and returns a pvalue for each gene')

parser.add_argument('--input', help='input h5ad', required=True)
parser.add_argument('--output1', help='output file for first bootstrapping step', required=True)
parser.add_argument('--output2', help='output file for second bootstrapping step', required=True)
parser.add_argument('--scenario', required=True, help='scenario of the data')

args = parser.parse_args()
file_in = args.input
file_out_bsc = args.output1
file_out_adv = args.output2
scenario = args.scenario

print(f'Input: {file_in}')
print(f'Output1: {file_out_bsc}')
print(f'Output2: {file_out_adv}')
print(f'Scenario: {scenario}')


def bootstrap(
    adata: AnnData,
    *,
    groupby: str,
    hierarchy: Sequence[str],
    sample_size: Tuple[int],
    aggr_func: Callable,
    use_raw=True,
    n=100,
) -> AnnData:
    """
    Perform hierarchical bootstrapping of single cells.

    Parameters
    ----------
    adata
        single cell data which may have hierarchical structure
    groupby
        The target column in `obs` with the groups between samples will be compared.
        E.g. leiden clusters, or treatment vs. control.
    hierarchy
        The levels of the hierarchical structure in the order from
        coarse to fine, e.g. `["dataset", "patient"]`.
    sample_size
        The number samples to draw from each level of the hierarchy.
        # TODO: The default is the smallest number of children a level of the
        hierarchy has.
        # TODO evaluate different strategies
    aggr_func
        Function used to aggregate the bootstrapping samples
    use_raw
        use data frame with raw data, stored in adata.raw
    n
        Number of bootstrap iterations

    Returns
    -------
    AnnData object with the same variables as the input, and `n` observation for
    each group. The generated anndata object can be used for comparisons without
    groups without having to worry about the bias introduced by pseudoreplication
    and the different levels of the hierarchical structure.
    """

    assert len(hierarchy) == len(sample_size)

    build_x = []
    build_obs = []

    X = adata.raw.X if use_raw else adata.X

    for i in tqdm(range(n), disable=False):
        for group in adata.obs[groupby].unique():
            group_mask = adata.obs[groupby] == group
            
            # Perform bootstrapping on level 1
            # Bootstrapping = sample with replacement, same size as original array
            # Following two lines are not incorporated in _get_idx() function since
            # we sample on a different size here
            unique_level1 = adata.obs[hierarchy[0]][group_mask].unique()
            level1_sample = np.random.choice(unique_level1, unique_level1.size)

            # Get index by recursively calling through all levels of the hierarchy
            idx = _get_idx(
                adata=adata,
                mask=group_mask,
                level_sample=level1_sample,
                hierarchy=hierarchy,
                sample_size=sample_size,
            )

            build_x.append(aggr_func(X[idx, :], axis=0))
            build_obs.append({"obs_names": f"{group}_{i}", groupby: group})

    return AnnData(
        X=np.vstack(build_x),
        var=adata.var,
        obs=pd.DataFrame().from_records(build_obs).set_index("obs_names"),
    )


def _get_idx(adata, mask, level_sample, hierarchy, sample_size):
    """
    Iterates over the levels of the hierarchy and draws the samples with replacement.

    Parameters
    ----------
    adata
        single cell data which may have hierarchical structure
    mask
        Mask of the last hierarchy level which will be refined in this method to the
        next level of the hierarchy
    level_sample
        contains drawn samples from the last hierarchy level. Method iterates over
        these and draws samples from the hierarchy level below.
    hierarchy
        The levels of the hierarchical structure in the order from
        coarse to fine, e.g. `["dataset", "patient"]`.
    sample_size
        The number samples to draw from each level of the hierarchy.
        TODO: The default is the smallest number of children a level of the
        hierarchy has.
    Returns
    -------
    Index of the drawn samples that can be applied to the adata object to generate
    bootstrapped adata
    """
    
    idx = []
    if len(hierarchy) > 1:
        for i in level_sample:
            level_mask = (adata.obs[hierarchy[0]] == i) & mask

            unique_next_level = adata.obs[hierarchy[1]][level_mask].unique()
            next_level_sample = np.random.choice(unique_next_level, sample_size[0])

            idx.extend(
                _get_idx(
                    adata=adata,
                    mask=level_mask,
                    level_sample=next_level_sample,
                    hierarchy=hierarchy[1:],
                    sample_size=sample_size[1:],
                )
            )

    else:
        for i in level_sample:
            level_mask = (adata.obs[hierarchy[0]] == i) & mask
            cells = np.where(level_mask)[0]
            idx.extend(np.random.choice(cells, sample_size[0]))

    return idx


def eval_hb(adata, var):
    X = adata.X

    mask_var = adata.var.index == var
    one_over_two = 0
    two_over_one = 0

    # Iterate over bootstrap samples
    for i in range(0, len(adata.obs_names), 2):
        if X[i, mask_var] > X[i + 1, mask_var]:
            one_over_two += 1
        else:
            two_over_one += 1

    return min(one_over_two / (adata.shape[0] / 2), two_over_one / (adata.shape[0] / 2))


adata = sc.read_h5ad(file_in)

sc.pp.normalize_total(adata, target_sum=1e6)

if scenario in ["atlas", "atlas_hvg", "atlas-less-de"]:
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Batch', 'Sample'],
                         sample_size=[10, 200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=10_000)
elif scenario in ["dataset", "dataset_hvg", "dataset-less-de"]:
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Sample'],
                         sample_size=[200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=10_000)
elif scenario in ["atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de"]:
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Batch', 'Sample'],
                         sample_size=[10, 200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=10_000)
elif scenario == "atlas-negative":
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Batch', 'Sample'],
                         sample_size=[10, 200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=10_000)
elif scenario == "kang2018":
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Batch', 'Sample'],
                         sample_size=[10, 200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=10_000)
elif scenario in ["dataset-ub-cells", "dataset-ub-cells_hvg", "dataset-ub-cells-less-de"]:
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Sample'],
                         sample_size=[200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=10_000)
elif scenario == "luca":
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Sample'],
                         sample_size=[200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=10_000)

p_values = {}

for gene in tqdm(adata_hb.var_names):
    p_values[gene] = eval_hb(adata_hb, gene)
p_values = pd.Series(p_values, name='pvalue')

p_values.to_csv(file_out_bsc, sep = '\t')
del p_values
del adata
del adata_hb


# Advanced bootstrapping starts here

pvalues = pd.read_csv(file_out_bsc, sep='\t', index_col=0)
pvalues_old = pvalues['pvalue']
test_genes = pvalues[pvalues_old == 0].index

adata = sc.read_h5ad(file_in)

sc.pp.normalize_total(adata, target_sum=1e6)

# Subset adata to necessary genes only
adata = adata[:,np.isin(adata.var.index, test_genes)]

if scenario in ["atlas", "atlas_hvg", "atlas-less-de"]:
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Batch', 'Sample'],
                         sample_size=[10, 200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=100_000)
elif scenario in ["dataset", "dataset_hvg", "dataset-less-de"]:
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Sample'],
                         sample_size=[200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=100_000)
elif scenario in ["atlas-ub-conditions", "atlas-ub-conditions_hvg", "atlas-ub-conditions-less-de"]:
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Batch', 'Sample'],
                         sample_size=[10, 200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=100_000)
elif scenario == "atlas-negative":
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Batch', 'Sample'],
                         sample_size=[10, 200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=100_000)
elif scenario == "kang2018":
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Batch', 'Sample'],
                         sample_size=[10, 200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=100_000)
elif scenario in ["dataset-ub-cells", "dataset-ub-cells_hvg", "dataset-ub-cells-less-de"]:
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Sample'],
                         sample_size=[200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=100_000)
elif scenario == "luca":
    adata_hb = bootstrap(adata=adata,
                         groupby='Condition',
                         hierarchy=['Sample'],
                         sample_size=[200],
                         aggr_func=np.mean,
                         use_raw=False,
                         n=100_000)
else:
    raise Exception("scenario not found")


pvalues_new = {}

for gene in tqdm(adata_hb.var_names):
    pvalues_new[gene] = eval_hb(adata_hb, gene)
pvalues_new = pd.Series(pvalues_new, name='pvalue_new')

pvalues['pvalue_new'] = pvalues_new

pvalues['pvalue'] = np.where(pvalues['pvalue'] == 0, pvalues['pvalue_new'], pvalues['pvalue'])
pvalues = pvalues.drop(columns=['pvalue_new'])

pvalues.to_csv(file_out_adv, sep = '\t')
