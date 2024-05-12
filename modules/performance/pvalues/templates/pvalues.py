#!/usr/bin/env python3

import pandas as pd
import numpy as np


meta_string = "${meta_string}"
path_string = "${path_string}"
file_out = "pval_${meta.scenario}_${meta.run}.tsv"


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
            res = pd.read_csv(path, sep='\\t', index_col=0)['pvalue']
        case 'distinct':
            res = pd.read_csv(path, sep='\\t', index_col=0)['p_val'].rename('pvalue')
        case 'mast':
            output = pd.read_csv(path, sep = '\\t')
            pval = output[np.logical_and(output['component'] == 'H',output['contrast'] == 'ConditionCondition2')]
            logfc = output[np.logical_and(output['component'] == 'logFC', output['contrast'] == 'ConditionCondition2')]
            res = pd.merge(pval, logfc, on='primerid', suffixes=('_pval', '_logfc')).set_index('primerid')
            res = res['Pr(>Chisq)_pval'].rename('pvalue')
        case 'scvi':
            res = pd.read_csv(path, sep = '\\t', index_col=0)['proba_not_de']
        case 'dream':
            res = pd.read_csv(path, sep = '\\t', index_col=0)['P.Value']
        case _:
            raise ValueError(f'Wrong method and file: {method} / {path}')

    pvalues[method] = res


results = pd.DataFrame.from_dict(pvalues)
results.to_csv(file_out, sep = "\\t")
