#!/usr/bin/env python3

import pandas as pd

meta = "${meta_string}".split(';')

# Process meta item into list of dictionaries
dicts = []
for item in meta:
    item = item.strip('[]')
    dictionary = {}
    for key_value_pair in item.split(', '):
        key, value = key_value_pair.split(':')
        dictionary[key] = value
    dicts.append(dictionary)

# Create dataframe from list of dictionaries
df = pd.DataFrame(dicts)
fixed_genes = df[df['scenario'] == 'fixed_genes'].drop(columns=['scenario', 'run']).sort_values(by=['n_cells', 'n_genes'])
fixed_cells = df[df['scenario'] == 'fixed_cells'].drop(columns=['scenario', 'run']).sort_values(by=['n_cells', 'n_genes'])

# Reorder columns of the dataframe
fixed_genes = fixed_genes[['n_cells', 'n_genes'] + fixed_genes.columns[2:].tolist()]
fixed_cells = fixed_cells[['n_cells', 'n_genes'] + fixed_cells.columns[2:].tolist()]

# Cast columns to numeric
assert fixed_genes.columns.tolist() == fixed_cells.columns.tolist()
for column in fixed_genes.columns[2:]:
    fixed_genes[column] = pd.to_numeric(fixed_genes[column], errors='raise')
    fixed_cells[column] = pd.to_numeric(fixed_cells[column], errors='raise')

# Compute mean on grouped dataframe
fixed_genes = fixed_genes.groupby(['n_cells', 'n_genes']).mean().reset_index()
fixed_cells = fixed_cells.groupby(['n_cells', 'n_genes']).mean().reset_index()

# Write dataframes to file
fixed_genes.to_csv('fixed_genes.tsv', sep='\\t', index=False)
fixed_cells.to_csv('fixed_cells.tsv', sep='\\t', index=False)