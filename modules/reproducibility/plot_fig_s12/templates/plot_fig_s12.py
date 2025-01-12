#!/usr/bin/env python3

import os
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from natsort import natsorted
from tqdm import tqdm


DATA_DIR = "."
N_FIXED = 100

METHOD_ORDER = [
    'deseq2',
    'dream',
    'hierarchical-bootstrapping',
    'mast',
    'permutation-test',
    'scvi',
    'distinct',
    'scdd',
    'ttest'
]

METHOD_TO_LEGEND = {
    'deseq2': 'DESeq2',
    'dream': 'DREAM',
    'hierarchical-bootstrapping': 'Hierarchical\\nBootstrapping',
    'mast': 'MAST',
    'permutation-test': 'Permutation\\nTest',
    'scvi': 'scVI',
    'distinct': 'distinct',
    'scdd': 'scDD',
    'ttest': 't-test'
}

# okabeito palette
PALETTE = [
    "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7", "#999999", "#000000"
]

files = natsorted([f for f in os.listdir(DATA_DIR) if f.endswith(".tsv")])
df_dict = {}
for f in files:
    df = pd.read_csv(os.path.join(DATA_DIR, f), sep="\t", index_col=0)
    df_dict[f] = df




all_genes = set(df_dict[files[0]].index)
for f in files[1:]:
    all_genes = all_genes.intersection(df_dict[f].index)

all_genes = sorted(all_genes)
for f in files:
    df_dict[f] = df_dict[f].loc[all_genes]



methods = df_dict[files[0]].columns
num_datasets = len(files)
num_genes = len(all_genes)


# Heatmap

combinations = [(method, f) for method in methods for f in files]
order_index = {method: index for index, method in enumerate(METHOD_ORDER)}
combinations = sorted(combinations, key=lambda x: order_index.get(x[0], len(METHOD_ORDER)))


all_top_sets = {}
for method, f in combinations:
    df = df_dict[f]
    top_genes = df[method].sort_values().head(N_FIXED).index
    all_top_sets[(method, f)] = set(top_genes)

n_combos = len(combinations)
jaccard_matrix = np.zeros((n_combos, n_combos))

for i in range(n_combos):
    for j in range(n_combos):
        set_i = all_top_sets[combinations[i]]
        set_j = all_top_sets[combinations[j]]
        intersection = set_i.intersection(set_j)
        union = set_i.union(set_j)
        jaccard_matrix[i, j] = len(intersection) / len(union) if union else 0


row_labels = [f"{m}-{os.path.splitext(f)[0].split('_')[2]}" for m, f in combinations]
df_jaccard = pd.DataFrame(jaccard_matrix, index=row_labels, columns=row_labels)


annotations_df = pd.DataFrame({
    "Method": [m for m, _ in combinations],
    "Dataset": [os.path.splitext(f)[0].split('_')[2] for _, f in combinations]
}, index=row_labels)

unique_methods = annotations_df["Method"].unique()
method_to_color = dict(zip(unique_methods, PALETTE))

row_colors = annotations_df["Method"].map(method_to_color)
col_colors = row_colors

g = sns.clustermap(
    df_jaccard,
    row_colors=row_colors,
    col_colors=col_colors,
    row_cluster=False,
    col_cluster=False,
    cmap="viridis",
    xticklabels=False,
    yticklabels=False,
    figsize=(20, 20),  # Increased figure width
    cbar_pos=(0.08, 0.1, 0.02, 0.6)  # [left, bottom, width, height]
)

# Customize colorbar
colorbar = g.ax_heatmap.collections[0].colorbar
colorbar.ax.tick_params(labelsize=16)

# Mapping labels to indices
index_mapping = {label: i for i, label in enumerate(row_labels)}
method_groups = annotations_df.groupby("Method").groups  # dict: method -> list of row_labels

# Adjust color axes
g.ax_row_colors.set_ylim(0, n_combos)
g.ax_row_colors.invert_yaxis()
g.ax_col_colors.set_xlim(0, n_combos)

# Remove ticks from color bars that group heatmap
g.ax_row_colors.set_xticks([])
g.ax_row_colors.set_yticks([])
g.ax_col_colors.set_xticks([])
g.ax_col_colors.set_yticks([])

# Remove whitespace reserved for dendrogram
g.ax_col_dendrogram.set_visible(False)

# Add method labels to color strips
for method, label_list in method_groups.items():
    # Convert row labels to integer positions:
    row_positions = [index_mapping[label] for label in label_list]
    min_r, max_r = min(row_positions), max(row_positions)
    center_r = (min_r + max_r + 1) / 2.0  # center of that block

    # Place text in row_colors strip:
    g.ax_row_colors.text(
        0.5,               # x ~ 0.5 to center horizontally in the color strip
        center_r,          # y is the vertical center
        METHOD_TO_LEGEND.get(method, method),
        ha="center",
        va="center",
        rotation=90,
        color="white",
        fontsize=12,
        fontweight="bold"
    )

    # For columns, same idea:
    min_c, max_c = min(row_positions), max(row_positions)
    center_c = (min_c + max_c + 1) / 2.0
    g.ax_col_colors.text(
        center_c,
        0.5,
        METHOD_TO_LEGEND.get(method, method),
        ha="center",
        va="center",
        rotation=0,
        color="white",
        fontsize=12,
        fontweight="bold"
    )

g.savefig("Fig_S12.png")
