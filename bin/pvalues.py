#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(
                    prog='Aggregate pvalues',
                    description='aggregates pvalues for every method and run')

parser.add_argument('--meta_string', help='string with the meta information from nextflow', required=True)
parser.add_argument('--path_string', help='string with the paths from nextflow', required=True)
parser.add_argument('--output', help='output file with pvalues [.tsv]', required=True)

args = parser.parse_args()
meta_string = args.meta_string
path_string = args.path_string
file_out = args.output


def meta_string_to_dict(meta: str) -> dict:
    meta = meta.strip('[]').replace(' ', '')
    meta_dict = {}
    for el in meta.split(','):
        key, value = el.split(':')
        meta_dict[key] = value
    return meta_dict


pvalues = {}
for meta, path in zip(meta_string.split(';'), path_string.split(';')):
    meta = meta_string_to_dict(meta)

    method = meta['method']
    match method:
        case 'hierarchical-bootstrapping' | 'permutation-test' | 'deseq2':
            res = pd.read_csv(path, sep='\t', index_col=0)['pvalue']
        case 'distinct':
            res = pd.read_csv(path, sep='\t', index_col=0)['p_val'].rename('pvalue')
        case 'mast':
            output = pd.read_csv(path, sep = '\t')
            pval = output[np.logical_and(output['component'] == 'H',output['contrast'] == 'ConditionCondition2')]
            logfc = output[np.logical_and(output['component'] == 'logFC', output['contrast'] == 'ConditionCondition2')]
            res = pd.merge(pval, logfc, on='primerid', suffixes=('_pval', '_logfc')).set_index('primerid')
            res = res['Pr(>Chisq)_pval'].rename('pvalue')
        case 'scvi':
            res = pd.read_csv(path, sep = '\t', index_col=0)['proba_not_de']
        case 'dream':
            res = pd.read_csv(path, sep = '\t', index_col=0)['P.Value']
        case _:
            raise ValueError(f'Wrong method and file: {method} / {path}')

    pvalues[method] = res


results = pd.DataFrame.from_dict(pvalues)
results.to_csv(file_out, sep = "\t")
