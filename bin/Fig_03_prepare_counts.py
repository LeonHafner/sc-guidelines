#!/usr/bin/env python3

import argparse
import scanpy as sc


parser = argparse.ArgumentParser(description="Prepare cells per sample counts")
parser.add_argument("--input", "--paths", required=True, help="String of file paths formatted as a python list")
parser.add_argument("--output", required=True, help="path of the output file")

args = parser.parse_args()
input = args.input
output = args.output


def parse_path_string(paths_string):
    """
    Parse the path_string into a list of strings.
    """
    return paths_string.strip('[]').split(', ')


def write_list_to_file(file_path, my_list):
    try:
        with open(file_path, 'w') as file:
            for item in my_list:
                file.write(str(item) + '\n')
        print("Elements written to the file successfully.")
    except IOError:
        print(f"Error writing to the file: {file_path}")


cells_per_sample = []
for file in parse_path_string(input):
    ad = sc.read(file)
    cells_per_sample += ad.obs.Sample.value_counts().tolist()


write_list_to_file(output, cells_per_sample)

