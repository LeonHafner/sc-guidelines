#!/usr/bin/env python3
import argparse
import scanpy as sc


parser = argparse.ArgumentParser(
                    prog='Preprocessing Step',
                    description='Program takes an adata file and removes all genes that are expressed in less than a specified amount of the cells')

parser.add_argument('--input', '-i', help='input h5ad', required=True)
parser.add_argument('--output', '-o', help='output h5ad', required=True)
parser.add_argument('--threshold', default=0.1, type=float, required=False, help='threshold of cells that have to express a gene in order to keep this gene in the data [0 - 1]')

args = parser.parse_args()
file_in = args.input
file_out = args.output
threshold = args.threshold

print(f'Input: {file_in}')
print(f'Output: {file_out}')
print(f'Threshold: {threshold}')

# Remove genes that are expressed in less than threshold of cells
adata = sc.read_h5ad(file_in)
adata = adata[:, (adata.X > 0).sum(axis=0) > (threshold * adata.n_obs)]

adata.write_h5ad(file_out)
