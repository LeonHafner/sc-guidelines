#!/usr/bin/env python
# coding: utf-8
import sys
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from sklearn.metrics import confusion_matrix, auc


parser = argparse.ArgumentParser(
                    prog='Precision-Recall computation',
                    description='Program computes the PRCs for a given tsv of pvalues (columns are methods)')

parser.add_argument('--pvalues', help='input pvalues [.tsv]', required=True)
parser.add_argument('--de_genes', help='de-genes for the dataset', required=True)
parser.add_argument('--output', help = "Output file", required=True)
parser.add_argument('--run', help = "Specify run", required=True)

args = parser.parse_args()
path_pvalues = args.pvalues
path_degenes = args.de_genes
output = args.output
run = args.run

print(f'p-values: {path_pvalues}')
print(f'de-genes: {path_degenes}')
print(f'Output: {output}')
print(f'Run: {run}')


pvalues = pd.read_csv(path_pvalues, sep='\t', index_col=0)
de_genes = pd.read_csv(path_degenes, sep='\t', index_col=0)


auc_methods = {'method': [], 'auc': []}
for method in pvalues:
    cut_offs = pvalues[method].sort_values().unique().tolist()

    precision = []
    recall = [0]

    y_true = np.isin(pvalues.index.tolist(), de_genes.index.tolist())
    for pvalue in tqdm(cut_offs[::10], disable = False):
        #y_pred = [1 if pvalues[method][gene] <= pvalue else 0 for gene in pvalues.index]
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

    pd.DataFrame.from_dict({'precision': precision, 'recall': recall}).to_csv(f'{output}/prc_{method}_{run}.tsv', sep='\t', index=False)

pd.DataFrame(auc_methods).to_csv(f'{output}/auc_{run}.tsv', sep = '\t', index=False)

