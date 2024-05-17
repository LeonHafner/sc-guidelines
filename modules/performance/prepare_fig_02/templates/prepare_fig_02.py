#!/usr/bin/env python3

import scanpy as sc


input = "${path_list}"

def parse_path_string(paths_string):
    """
    Parse the path_string into a list of strings.
    """
    return paths_string.split(' ')


def write_list_to_file(file_path, my_list):
    try:
        with open(file_path, 'w') as file:
            for item in my_list:
                file.write(str(item) + '\\n')
        print('Elements written to the file successfully.')
    except IOError:
        print(f'Error writing to the file: {file_path}')


cells_per_sample = []
for file in parse_path_string(input):
    ad = sc.read(file)
    cells_per_sample += ad.obs.Sample.value_counts().tolist()


write_list_to_file('cell_counts.txt', cells_per_sample)

