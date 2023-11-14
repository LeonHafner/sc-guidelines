#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
from os.path import join, exists, basename
import argparse

parser = argparse.ArgumentParser(
                    prog='Aggregate pvalues',
                    description='aggregates pvalues for every method and run')

parser.add_argument('--input', help='directory containing a directory for each method', required=True)
parser.add_argument('--output', help='output file with pvalues [.tsv]', required=True)
parser.add_argument('--run', help='specify run to process', required=False)

args = parser.parse_args()
dir_in = args.input
file_out = args.output
run = args.run


print(f'Input: {dir_in}')
print(f'Output: {file_out}')
print(f'Run: {run}')


poss_files = [
    f'5_mast/mast_{run}.tsv',
    f'6_distinct/distinct_{run}.tsv',
    f'7_deseq2/deseq2_{run}.tsv',
    f'8_permutation/permutation_{run}.tsv',
    f'9_hierarchical-bootstrapping/hierarchical-bootstrapping_adv_{run}.tsv',
    f'10_scvi/scvi_{run}.tsv',
    f'11_dream/dream_{run}.tsv'
]


files = []
for file in poss_files:
    if exists(join(dir_in, file)):
        files.append(join(dir_in, file))
    else:
        print(f'{file} missing')


pvalues = {}

for file in files:
    method = '_'.join(basename(file).split('.')[0].split('_')[:-1])
    if method == 'hierarchical-bootstrapping_adv':
        method = 'hierarchical-bootstrapping'
        res = pd.read_csv(file, sep='\t', index_col=0)['pvalue']
    elif method == 'permutation':
        res = pd.read_csv(file, sep='\t', index_col=0)['pvalue']
    elif method == 'deseq2':
        res = pd.read_csv(file, sep='\t', index_col=0)['pvalue']
    elif method == 'distinct':
        res = pd.read_csv(file, sep='\t', index_col=0)['p_val'].rename('pvalue')
        res.index.name = None
    elif method == 'mast':
        output = pd.read_csv(file, sep = '\t')
        pval = output[np.logical_and(output['component'] == 'H',output['contrast'] == 'ConditionCondition2')]
        logfc = output[np.logical_and(output['component'] == 'logFC', output['contrast'] == 'ConditionCondition2')]
        res = pd.merge(pval, logfc, on='primerid', suffixes=('_pval', '_logfc')).set_index('primerid')
        res = res['Pr(>Chisq)_pval'].rename('pvalue')
    elif method == 'scvi':
        res = pd.read_csv(file, sep = '\t', index_col=0)['proba_not_de']
    elif method == 'dream':
        res = pd.read_csv(file, sep = '\t', index_col=0)['P.Value']
    else:
        raise Exception('Wrong file', file)

    pvalues[method] = res


results = pd.DataFrame.from_dict(pvalues)
results.to_csv(file_out, sep = "\t")

