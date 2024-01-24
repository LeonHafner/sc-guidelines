#!/usr/bin/env python3

import scanpy as sc
import numpy as np
import argparse
import random



parser = argparse.ArgumentParser()

parser.add_argument('--input', '-i', help='input [.h5ad]', required=True)
parser.add_argument('--output', '-o', help='output [.h5ad]', required=True)

args = parser.parse_args()
file_in = args.input
file_out = args.output


print(f'Input: {file_in}')
print(f'Output: {file_out}')


ad = sc.read(file_in)


# Subset samples of specified batch and condition
samples_batch1_condition1 = ad.obs[np.logical_and(ad.obs["Batch"] == "Batch1", ad.obs["Condition"] == "Condition1")]["Sample"].unique().tolist()
samples_batch1_condition2 = ad.obs[np.logical_and(ad.obs["Batch"] == "Batch1", ad.obs["Condition"] == "Condition2")]["Sample"].unique().tolist()
samples_batch2_condition1 = ad.obs[np.logical_and(ad.obs["Batch"] == "Batch2", ad.obs["Condition"] == "Condition1")]["Sample"].unique().tolist()
samples_batch2_condition2 = ad.obs[np.logical_and(ad.obs["Batch"] == "Batch2", ad.obs["Condition"] == "Condition2")]["Sample"].unique().tolist()


# Draw from samples with specified batch and condition
drawn_samples_batch1_condition1 = random.sample(samples_batch1_condition1, k=9)
drawn_samples_batch1_condition2 = random.sample(samples_batch1_condition2, k=1)
drawn_samples_batch2_condition1 = random.sample(samples_batch2_condition1, k=1)
drawn_samples_batch2_condition2 = random.sample(samples_batch2_condition2, k=9)

# Subset anndata to drawn samples
drawn_samples = drawn_samples_batch1_condition1 + drawn_samples_batch1_condition2 + drawn_samples_batch2_condition1 + drawn_samples_batch2_condition2
ad = ad[np.isin(ad.obs['Sample'], drawn_samples)]

ad.write_h5ad(file_out)