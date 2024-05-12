#!/usr/bin/env python3

import pandas as pd
import numpy as np
from tqdm import tqdm

pvalues = pd.read_csv("${pvalues}", sep='\\t', index_col=0)
scvi_values = pd.read_csv("${scvi_values}", sep='\\t')

# Merge pvals from scvi with other methods
cut_offs = sorted(list(set(pvalues.stack().dropna().tolist() + scvi_values['pvals'].tolist())))

results_fpr = {'pvals': cut_offs}
results_fps = {'pvals': cut_offs}
# For all methods except scvi
for method in pvalues:
    fpr = []
    fps = []
    
    for pvalue in tqdm(cut_offs):
        fp = np.sum(pvalues[method].to_numpy() <= pvalue)
        tn = len(pvalues[method]) - fp
        fpr.append(fp / (fp + tn))
        fps.append(fp)

    results_fpr[method] = fpr
    results_fps[method] = fps

# For scvi
fpr = []
fps = []
for pvalue in cut_offs:
    fp = scvi_values[scvi_values['pvals'] <= pvalue].sort_values(by='pvals', ascending=False).iloc[0, 0]
    tn = scvi_values['fps'].max() - fp
    fpr.append(fp / (fp + tn))
    fps.append(fp)
results_fpr['scvi'] = fpr
results_fps['scvi'] = fps

pd.DataFrame(results_fpr).to_csv("${meta.scenario}_${meta.run}_fpr.tsv", sep = "\\t", index = False)
pd.DataFrame(results_fps).to_csv("${meta.scenario}_${meta.run}_fps.tsv", sep = "\\t", index = False)