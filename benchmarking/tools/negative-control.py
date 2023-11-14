#!/usr/bin/env python
# coding: utf-8
import sys
import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from sklearn.metrics import confusion_matrix

parser = argparse.ArgumentParser(
                    prog='false-positive count computation')

parser.add_argument('--pvalues_adj', help='input pvalues adjusted [.tsv]', required=True)
parser.add_argument('--output', help = "Output file", required=True)

args = parser.parse_args()
path_pvalues = args.pvalues
output = args.output


print(f'p-values: {path_pvalues}')
print(f'Output: {output}')


pvalues = pd.read_csv(path_pvalues, sep='\t', index_col=0)
cut_offs = sorted(list(set(pvalues.stack().dropna().tolist())))


prc = {}
for method in pvalues:
    print(method)
    false_positives = []
    
    for pvalue in tqdm(cut_offs):
        false_positives.append(np.sum(pvalues[method].to_numpy() <= pvalue))
        
    pd.DataFrame.from_dict({'cutoff': cut_offs, 'false_positives': false_positives}).to_csv(f'{output}/fps_{method}.tsv', sep='\t', index=False)

