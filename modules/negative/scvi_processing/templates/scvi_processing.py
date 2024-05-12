#!/usr/bin/env python3

import os
import pandas as pd

entries = [os.path.join('result_files', entry) for entry in os.listdir('result_files') if entry.startswith('result-')]

list_of_series = []
for entry in entries:
    result_table = pd.read_csv(entry, sep="\\t", index_col=0)
    # iloc transforms from pd.DataFrame into pd.Series
    result_series = result_table.loc[:,result_table.columns.str.startswith("is_de_fdr")].iloc[:, 0]
    list_of_series.append(result_series)

df = pd.concat(list_of_series, axis=1).sort_index(axis=1)
df = pd.DataFrame({'fps': df.sum(axis=0)})
df['pvals'] = [float(element.split('_')[3]) for element in df.index.tolist()]
df.to_csv("scvi_values.tsv", sep="\\t", index=False)