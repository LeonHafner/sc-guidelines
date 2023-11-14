#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import numpy as np
import argparse
import random

parser = argparse.ArgumentParser()

parser.add_argument('--input', help='input [.h5ad]', required=True)
parser.add_argument('--output', help='output [.h5ad]', required=True)

args = parser.parse_args()
file_in = args.input
file_out = args.output

print(f'Input: {file_in}')
print(f'Output: {file_out}')



ad = sc.read(file_in)

# Exclude 6 Samples from Batch1 and Condition1
excl1 = random.sample(ad.obs[np.logical_and(ad.obs['Batch'] == 'Batch1', ad.obs['Condition'] == "Condition1")]['Sample'].drop_duplicates().tolist(), 6)

# Exclude 6 Samples from Batch2 and Condition2
excl2 = random.sample(ad.obs[np.logical_and(ad.obs['Batch'] == 'Batch2', ad.obs['Condition'] == "Condition2")]['Sample'].drop_duplicates().tolist(), 6)

ad = ad[np.logical_not(np.isin(ad.obs['Sample'], excl1 + excl2))]

ad.write_h5ad(file_out)
