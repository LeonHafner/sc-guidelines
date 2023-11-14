#!/usr/bin/env python3
import sys
import argparse
import scanpy as sc
import numpy as np

parser = argparse.ArgumentParser(
                    prog='DE-Genes acquiring',
                    description='Program extracts differentially expressed genes from the splimp ground truth')

parser.add_argument('--input', help='input h5ad', required=True)
parser.add_argument('--output', help='output tsv', required=True)

args = parser.parse_args()
file_in = args.input
file_out = args.output

print(f'Input: {file_in}')
print(f'Output: {file_out}')


adata = sc.read_h5ad(file_in)
de_genes = np.logical_or(adata.var['ConditionDE.Condition1'] != 1, adata.var['ConditionDE.Condition2'] != 1)
adata.var[['ConditionDE.Condition1', 'ConditionDE.Condition2']][de_genes].to_csv(file_out, sep = '\t')
