#!/usr/bin/env python3

import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.metrics import confusion_matrix, auc


path_pvalues = "${pvalues}"
path_ground_truth = "${ground_truth}"
output = "${meta.scenario}_${meta.run}"

print(f'p-values: {path_pvalues}')
print(f'ground truth: {path_ground_truth}')
print(f'Output: {output}')

# Load pvalues and ground truth
pvalues = pd.read_csv(path_pvalues, sep='\\t', index_col=0)
ground_truth = pd.read_csv(path_ground_truth, sep='\\t', index_col=0)

# Extract DE genes from the ground truth DataFrame
de_genes = ground_truth[np.logical_or(ground_truth['ConditionDE.Condition1'] != 1, ground_truth['ConditionDE.Condition2'] != 1)]

# Get genes that were removed by filtering before methods were applied, set their pvalue to 1 for all methods
removed_genes = ground_truth[np.logical_not(np.isin(ground_truth.index, pvalues.index))].index
for gene in removed_genes:
    pvalues.loc[gene] = [1.0] * len(pvalues.columns)

# Some methods predict pvalues NA
pvalues = pvalues.fillna(1.0)

auc_methods = {'method': [], 'auc': []}
for method in pvalues:
    cut_offs = pvalues[method].sort_values().unique().tolist()

    precision = []
    recall = [0]

    y_true = np.isin(pvalues.index.tolist(), de_genes.index.tolist())
    for pvalue in tqdm(cut_offs[::10], disable = False):
        y_pred = (pvalues[method].to_numpy() <= pvalue).astype(int)
        (TN, FP), (FN, TP) = confusion_matrix(y_true=y_true, y_pred=y_pred)
        precision.append(TP / (TP + FP))
        recall.append(TP / (TP + FN))

    precision = [precision[0]] + precision + [sum(y_true) / len(y_true)]
    recall.append(1)
    
    prc = pd.DataFrame.from_dict({'precision': precision, 'recall': recall})
    prc.dropna(inplace = True)
        
    area_under_curve = auc(prc.recall, prc.precision)
    auc_methods['method'].append(method)
    auc_methods['auc'].append(area_under_curve)

    pd.DataFrame.from_dict({'precision': precision, 'recall': recall}).to_csv(f'prc_{output}_{method}.tsv', sep='\\t', index=False)

pd.DataFrame(auc_methods).to_csv(f'auc_{output}.tsv', sep = '\\t', index=False)
