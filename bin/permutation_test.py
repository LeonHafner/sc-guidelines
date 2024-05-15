#!/usr/bin/env python3

import argparse
import scanpy as sc
import pandas as pd
from anndata import AnnData
from tqdm import tqdm
import numpy as np

np.random.seed(0)

parser = argparse.ArgumentParser(
                    prog='Permutation Approach',
                    description='Program works and pseudobulked data and evaluates pvalues by permuting group labels')

parser.add_argument('--input', help='input h5ad', required=True)
parser.add_argument('--output', help='output tsv', required=True)
parser.add_argument('-n', required=False, default=10_000, type=int, help = "Number of minimal iterations to perform")
parser.add_argument('-n_max', required=False, default=100_000, type=int, help = "Number of maximal iterations to perform")
parser.add_argument('--scenario', required=True, help='Scenario of the data')
    

args = parser.parse_args()
file_in = args.input
file_out = args.output
n = args.n
n_max = args.n_max
scenario = args.scenario

print(f'Input: {file_in}')
print(f'Output: {file_out}')
print(f'min iterations: {n}')
print(f'max iterations: {n_max}')


def permutation(
    adata: AnnData,
    groupby: str,
    group1: str,
    group2: str,
    var: str,
    n: int,
    n_max: int,
    use_raw=False,
) -> float:
    """
    Compute a p-value for two groups in a permutation based manner.
    Currently, only working for two differing groups.

    Parameters
    ----------
    adata
        single cell data
    groupby
        The target column in 'obs' with the groups between samples will be compared.
    group1
        The first group to be tested for significance.
    group2
        The second group to be tested for significance.
    var
        Variable examined for statistical significance.
    n
        number of iterations the group labels are permuted.
    n_max
        maximal number of iterations, only reached as long as p-value is still 0
    use_raw
        use data frame with raw data, stored in adata.raw

    Returns
    -------
    P-value, which indicates whether the difference in mean between the groups (
    'group1' and 'group2') in the variable 'var' is significant.
    The minimal achievable p-value is 1/n, which means that a vast amount of iterations
    is needed.
    """

    X = adata.raw.X if use_raw else adata.X

    mask_group1 = adata.obs[groupby] == group1
    mask_group2 = adata.obs[groupby] == group2
    
    mask_var = adata.var.index == var

    mean_group1 = np.mean(X[mask_group1, mask_var])
    mean_group2 = np.mean(X[mask_group2, mask_var])

    init_mean_diff = abs(mean_group1 - mean_group2)
    perm_mean_diffs = []
    perm_mean_greater = [0]
    iterations = 0

    while (perm_mean_greater is None) or (sum(perm_mean_greater) == 0 and iterations < n_max):
        iterations += n

        for _ in tqdm(range(n), disable=True):
            # Permute group labels
            perm_groups = np.random.permutation(adata.obs[groupby])

            # Mask group labels
            perm_mask_group1 = perm_groups == group1
            perm_mask_group2 = perm_groups == group2

            # Calculate the mean for each group
            mean_group1 = np.mean(X[perm_mask_group1, mask_var])
            mean_group2 = np.mean(X[perm_mask_group2, mask_var])
            perm_mean_diffs.append(abs(mean_group1 - mean_group2))

        # Count how often the difference in means is larger -> p-value
        perm_mean_greater = perm_mean_diffs > init_mean_diff
    return sum(perm_mean_greater) / len(perm_mean_greater)



adata = sc.read_h5ad(file_in)

group1 = "Condition1"
group2 = "Condition2"


p_values = {}
for gene in tqdm(adata.var.index):
    p_values[gene] = permutation(adata,
                                 groupby='Condition',
                                 group1=group1,
                                 group2=group2,
                                 var=gene,
                                 n=n,
                                 n_max=n_max)
p_values = pd.Series(p_values, name='pvalue')

p_values.to_csv(file_out, sep='\t')
